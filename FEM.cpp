#include <iostream>
#include <stdexcept>
#include <array>
#include <cmath>
#include "Headers.h"
#include "FEMAxisymmetric.h"
#include "FEMCartesian.h"

namespace {

	constexpr double kAxisymmetricTwoPi = 6.28318530717958647692;
	constexpr double kAxisymmetricRadiusEps = 1e-14;

	struct TriangleQuadraturePoint {
		double weight;
		std::array<double, 3> phi;
	};

	const std::array<TriangleQuadraturePoint, 3> kTriangleQuadrature = { {
		{ 1.0 / 3.0, { 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0 } },
		{ 1.0 / 3.0, { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 } },
		{ 1.0 / 3.0, { 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0 } }
	} };

	double signed_double_area(const Point& p1, const Point& p2, const Point& p3) {
		return (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y);
	}

}

std::ostream& operator<<(std::ostream& output, const Point& p) {
	output << "{ " << p.x << "; " << p.y << " }";
	return output;
}

Point operator+(const Point& a, const Point& b) {
	return { a.x + b.x, a.y + b.y };
}

Point operator-(const Point& a, const Point& b) {
	return { a.x - b.x, a.y - b.y };
}

Point operator*(const double a, const Point& p) {
	return { p.x * a, p.y * a };
}

Point operator/(const Point& p, const double a) {
	return { p.x / a, p.y / a };
}

function mul(function a, function b) {
	return [a, b](const Point& p) { return b(p) * a(p); };
}

function get_x(vec_function a) {
	return [a](const Point& p) { return a(p).x; };
}

function get_y(vec_function a) {
	return [a](const Point& p) { return a(p).y; };
}


FEM::FEM(const Point& a, const Point& b, size_t n, size_t m,
	CoordinateSystem coordinate_system)
	: mx(m), ny(n), A(2 * n * m, 2 * n * m), F(2 * n * m),
	left_down(a), right_up(b), coordinate_system(coordinate_system) {

	double dx = (b.x - a.x) / (m - 1), dy = (b.y - a.y) / (n - 1);
	u = std::vector<Point>(m * n, { NAN, NAN });
	triangle_area = dx * dy / 2;
	
	points.resize(m * n);
	f.resize(m * n);
	for (size_t i = 0; i != points.size(); ++i) {
		points[i] = { a.x + (i % m) * dx, a.y + (i / m) * dy };
	}

	triangles.resize(2 * (m - 1) * (n - 1));
	for (size_t i = 0; i != triangles.size(); ++i) {
		size_t x = i % (2 * (m - 1)), y = i / (2 * (m - 1));
		if (x % 2 == 0) { // "верхний" треугольник
			triangles[i] = {
				(y * m) + x / 2,
				((y + 1) * m) + x / 2,
				((y + 1) * m) + x / 2 + 1
			};
		}
		else { // "нижний" треугольник
			triangles[i] = {
				(y * m) + x / 2,
				((y + 1) * m) + x / 2 + 1,
				(y * m) + x / 2 + 1
			};
		}
	}
}

Point& FEM::get_point(size_t n) {
	return points[n];
}

const Point& FEM::get_point(size_t n) const {
	return points[n];
}

Triangle& FEM::get_triangle(size_t n) {
	return triangles[n];
}

const Triangle& FEM::get_triangle(size_t n) const {
	return triangles[n];
}

void FEM::print_points() const {
	for (size_t i = 0; i < points.size(); ++i) {
		std::cout << (*this)[i] << "\t";
	}
	std::cout << '\n';
}

void FEM::print_triangles() const {
	for (size_t i = 0; i < triangles.size(); ++i) {
		std::cout << "[ " <<
			(*this)(i).a << ", " <<
			(*this)(i).b << ", " <<
			(*this)(i).c << "]\n";
	}
}

void FEM::set_boundaries(char side, const vec_function& g) {
	if (side == 'S') { // Нижняя граница
		for (size_t j = 0; j != mx; ++j) {
			Point t = g(points[j]);
			if(isnan(u[j].x))
				u[j].x = t.x;
			if (isnan(u[j].y))
				u[j].y = t.y;
		}
	}
	else if (side == 'E') { // Правая граница
		for (size_t i = 0; i != ny; ++i) {
			Point t = g(points[mx - 1 + i * mx]);
			if (isnan(u[mx - 1 + i * mx].x))
				u[mx - 1 + i * mx].x = t.x;
			if (isnan(u[mx - 1 + i * mx].y))
				u[mx - 1 + i * mx].y = t.y;
		}
	}
	else if (side == 'N') { // Верхняя граница
		for (size_t j = 0; j != mx; ++j) {
			Point t = g(points[mx * (ny - 1) + j]);
			if (isnan(u[mx * (ny - 1) + j].x))
				u[mx * (ny - 1) + j].x = t.x;
			if (isnan(u[mx * (ny - 1) + j].y))
				u[mx * (ny - 1) + j].y = t.y;
		}
	}
	else { // Левая граница
		for (size_t i = 0; i != ny; ++i) {
			Point t = g(points[i * mx]);
			if (isnan(u[i * mx].x))
				u[i * mx].x = t.x;
			if (isnan(u[i * mx].y))
				u[i * mx].y = t.y;
		}
	}
}

void FEM::construct_AF(double E, double nu, vec_function f) {
	if (coordinate_system == CoordinateSystem::Cartesian) {
		fem_cartesian::assemble(A, F, points, triangles, triangle_area, E, nu, f);
	}
	else {
		fem_axisymmetric::assemble(A, F, points, triangles, triangle_area, E, nu, f);
	}
}

void FEM::apply_boundaries() {
	for (size_t i = 0; i < psize(); ++i) {
		if (!isnan(u[i].x)) {
			F[2 * i] = u[i].x;
			for (size_t j = 0; j < A[2 * i].size(); ++j)
				A[2 * i][j] = 0;
			A[2 * i][2 * i] = 1;
		}
		if (!isnan(u[i].y)) {
			F[2 * i + 1] = u[i].y;
			for (size_t j = 0; j < A[2 * i + 1].size(); ++j)
				A[2 * i + 1][j] = 0;
			A[2 * i + 1][2 * i + 1] = 1;
		}
	}
}

std::vector<Point> FEM::solve() {
	auto [L, U] = LU_decomposition(A);
	std::vector<double> u = solveLU(L, U, F);

	std::vector<Point> res(psize());
	for (size_t i = 0; i < psize(); ++i) {
		res[i].x = u[2 * i];
		res[i].y = u[2 * i + 1];
	}

	clear_AFu();

	return res;
}

void FEM::clear_AFu() {
	for (size_t i = 0; i < 2 * psize(); ++i) {
		F[i] = 0;
		if (i % 2 == 0)
			u[i / 2] = { NAN, NAN };
		for (size_t j = 0; j < 2 * psize(); ++j)
			A[i][j] = 0;
	}
}

std::pair<Matrix, std::vector<double>> FEM::get_AF() {
	return { A, F };
}

void FEM::set_AF(const Matrix& A_new, const std::vector<double>& F_new) {
	A = A_new;
	F = F_new;
}

void FEM::bc2_side(lambda_func j, size_t start, size_t finish, double len,
	int side, bool prev_side, bool next_side, 
	const std::vector<vec_function>& g, std::vector<double>& p_vec) {
	for (size_t i = start; i < finish; ++i) { // индекс по границе как в мксэ
		const Point midpoint = (points[j(i)] + points[j(i + 1)]) / 2;
		const double weight = (coordinate_system == CoordinateSystem::Cartesian)
			? fem_cartesian::boundary_segment_weight(midpoint, len)
			: fem_axisymmetric::boundary_segment_weight(midpoint, len);
		Point integral = weight * g[side](midpoint);

		// текущий узел
		p_vec[2 * j(i)] += integral.x;
		p_vec[2 * j(i) + 1] += integral.y;

		p_vec[2 * j(i + 1)] += integral.x;
		p_vec[2 * j(i + 1) + 1] += integral.y;
		
		/*if (!prev_side)
			p_vec[2 * j(start)] = p_vec[2 * j(start) + 1] = 0;

		if (!next_side)
			p_vec[2 * j(finish)] = p_vec[2 * j(finish) + 1] = 0;*/
	}
}

void FEM::calculate_bc2(const std::vector<size_t>& pos,
	const std::vector<vec_function>& g, std::vector<double>& p_vec) {
	
	double len_vert = (right_up.y - left_down.y) / (ny - 1),
		len_hor = (right_up.x - left_down.x) / (mx - 1);

	if (pos[0])  // слева ГУ 2 рода
		bc2_side([&](size_t i) { return i * mx; }, 0, ny - 1, 
			len_vert, 0, pos[3], pos[1], g, p_vec);

	if (pos[1])  // сверху ГУ 2 рода
		bc2_side([&](size_t i) { return (mx - 1) * (ny - 1) + i; }, 
			ny - 1, ny + mx - 2, len_hor, 1, pos[0], pos[2], g, p_vec);

	if (pos[2])  // справа ГУ 2 рода
		bc2_side([&](size_t i) { return mx * (2 * ny + mx - i - 2) - 1; }, 
			ny + mx - 2, 2 * ny + mx - 3, len_vert, 
			2, pos[1], pos[3], g, p_vec);
	
	if (pos[3])  // снизу ГУ 2 рода
		bc2_side([&](size_t i) { return 1 + 2 * (ny + mx) - 5 - i; }, 
			2 * ny + mx - 3, 2 * (ny + mx) - 4, len_hor, 
			3, pos[2], pos[0], g, p_vec);

}
