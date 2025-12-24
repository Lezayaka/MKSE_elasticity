#include <iostream>
#include <stdexcept>
#include "Headers.h"

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


FEM::FEM(const Point& a, const Point& b, size_t n, size_t m)
	: mx(m), ny(n), A(2 * n * m, 2 * n * m), F(2 * n * m),
	left_down(a), right_up(b){ // почему так m и n ?????????

	double dx = (b.x - a.x) / (m - 1), dy = (b.y - a.y) / (n - 1);
	u = std::vector<Point>(m * n, { NAN, NAN });
	triangle_area = dx * dy / 2;
	
	points.resize(m * n);
	f.resize(m * n);
	for (size_t i = 0; i != points.size(); ++i) {
		points[i] = { (i % m) * dx, (i / m) * dy };
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
	double lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
	double mu = E / (2 * (1 + nu));
	Matrix C({
		{lambda + 2 * mu, lambda, 0},
		{lambda, lambda + 2 * mu, 0},
		{0, 0, mu}
	});

	for (size_t t = 0; t != tsize(); ++t) {
		Triangle& T = triangles[t];
		Point &p1 = points[T.a], &p2 = points[T.b], &p3 = points[T.c];
		Matrix grad_phi({
			{(p2.y - p3.y) / (2 * triangle_area), (p3.x - p2.x) / (2 * triangle_area)},
			{(p3.y - p1.y) / (2 * triangle_area), (p1.x - p3.x) / (2 * triangle_area)},
			{(p1.y - p2.y) / (2 * triangle_area), (p2.x - p1.x) / (2 * triangle_area)}
		});

		Matrix R(3, 6ull);
		for (size_t p = 0; p < 3; ++p) {
			R[0][2 * p] = grad_phi[p][0];
			R[1][2 * p + 1] = grad_phi[p][1];

			R[2][2 * p] = grad_phi[p][1];
			R[2][2 * p + 1] = grad_phi[p][0];
		}
		
		Matrix Ae = triangle_area * R.T().dot(C).dot(R);
		

		std::vector<double> Fe(6, 0);
		Point center = (p1 + p2 + p3) / 3;

		Fe[0] = get_x(f)(center) * triangle_area / 3;
		Fe[1] = get_y(f)(center) * triangle_area / 3;

		Fe[2] = get_x(f)(center) * triangle_area / 3;
		Fe[3] = get_y(f)(center) * triangle_area / 3;

		Fe[4] = get_x(f)(center) * triangle_area / 3;
		Fe[5] = get_y(f)(center) * triangle_area / 3;
		

		auto p = [T](size_t i) {
			if (i == 0) return 2 * T.a;
			if (i == 1) return 2 * T.a + 1;
			if (i == 2) return 2 * T.b;
			if (i == 3) return 2 * T.b + 1;
			if (i == 4) return 2 * T.c;
			if (i == 5) return 2 * T.c + 1;
		};

		for (size_t i = 0; i < 6; ++i) {
			F[p(i)] += Fe[i];
			for (size_t j = 0; j < 6; ++j)
				A[p(i)][p(j)] += Ae[i][j];
		}

	}
	
	//A.print();
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

	clear_AF();

	return res;
}

void FEM::clear_AF() {
	for (size_t i = 0; i < 2 * psize(); ++i) {
		F[i] = 0;
		if (i % 2 == 0)
			u[i / 2] = { NAN, NAN };
		for (size_t j = 0; j < 2 * psize(); ++j)
			A[i][j] = 0;
	}
}

void FEM::bc2_side(lambda_func j, int start, int finish, double len,
	int side, bool prev_side, bool next_side, 
	const std::vector<vec_function>& g, std::vector<double>& p_vec) {
	for (int i = start; i < finish; ++i) { // индекс по границе как в мксэ

		Point integral = 0.5 * len *
			g[side]((points[j(i)] + points[j(i + 1)]) / 2);

		// текущий узел
		p_vec[2 * j(i)] += integral.x;
		p_vec[2 * j(i) + 1] += integral.y;

		p_vec[2 * j(i + 1)] += integral.x;
		p_vec[2 * j(i + 1) + 1] += integral.y;
		
		if (!prev_side)
			p_vec[2 * j(start)] = p_vec[2 * j(start) + 1] = 0;

		if (!next_side)
			p_vec[2 * j(finish)] = p_vec[2 * j(finish) + 1] = 0;
	}
}

void FEM::calculate_bc2(const std::vector<size_t>& pos,
	const std::vector<vec_function>& g, std::vector<double>& p_vec) {
	
	double len_vert = (right_up.y - left_down.y) / (ny - 1),
		len_hor = (right_up.x - left_down.x) / (mx - 1);

	if (pos[0])  // слева ГУ 2 рода
		bc2_side([&](int i) { return i * mx; }, 0, ny - 1, 
			len_vert, 0, pos[3], pos[1], g, p_vec);

	if (pos[1])  // сверху ГУ 2 рода
		bc2_side([&](int i) { return (mx - 1) * (ny - 1) + i; }, 
			ny - 1, ny + mx - 2, len_hor, 1, pos[0], pos[2], g, p_vec);

	if (pos[2])  // справа ГУ 2 рода
		bc2_side([&](int i) { return mx * (2 * ny + mx - i - 2) - 1; }, 
			ny + mx - 2, 2 * ny + mx - 3, len_vert, 
			2, pos[1], pos[3], g, p_vec);
	
	if (pos[3])  // снизу ГУ 2 рода
		bc2_side([&](int i) { return 1 + 2 * (ny + mx) - 5 - i; }, 
			2 * ny + mx - 3, 2 * (ny + mx) - 4, len_hor, 
			3, pos[2], pos[0], g, p_vec);

}