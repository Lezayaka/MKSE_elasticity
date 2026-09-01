#pragma once

// Core data structures and solvers for the finite superelement model.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

struct Point {
	double x, y;

	friend std::ostream& operator<<(std::ostream& output, const Point& p);
};

Point operator+(const Point& a, const Point& b);
Point operator-(const Point& a, const Point& b);
Point operator*(double a, const Point& p);
Point operator/(const Point& p, double a);

std::vector<double> operator-(const std::vector<double>& a,
	const std::vector<double>& b);

using function = std::function<double(const Point&)>;
using vec_function = std::function<Point(const Point&)>;
using lambda_func = std::function<size_t(size_t)>;

enum class CoordinateSystem {
	Cartesian,
	Axisymmetric
};

enum class ContactMethod {
	SlaveNodes,
	UniformLambdaPartition,
	UniformUnionPartition
};

enum class ContactSlaveBody {
	Bottom,
	Top
};

struct ContactOptions {
	CoordinateSystem coordinate_system = CoordinateSystem::Axisymmetric;
	ContactMethod method = ContactMethod::UniformUnionPartition;
	size_t lambda_node_count = 0;
	ContactSlaveBody slave_body = ContactSlaveBody::Bottom;
};

std::string to_string(CoordinateSystem coordinate_system);
std::string to_string(ContactMethod contact_method);
std::string to_string(ContactSlaveBody slave_body);

struct Triangle {
	size_t a, b, c;
};

struct Rectangle {
	size_t a, b, c, d;
};

struct MortarElement {
	double chi_left;
	double chi_right;
};

struct ContactDiscretization {
	std::vector<double> lambda_nodes;
	std::vector<MortarElement> mortar_elements;
};

class Matrix {
private:
	std::vector<std::vector<double>> matrix;

public:
	Matrix(size_t n, double a = 0);
	Matrix(size_t n, size_t m, double a = 0);
	Matrix(std::vector<std::vector<double>> m) : matrix(m) {};
	Matrix(Matrix* M) { matrix = M->matrix; };

	std::vector<double>& operator[](size_t i) { return matrix[i]; };
	const std::vector<double>& operator[](size_t i) const { return matrix[i]; };
	Matrix operator*(double a);
	Matrix& operator*=(double a);

	Matrix T() const;
	size_t size(short axis = 0) const {
		return axis ? matrix[0].size() : matrix.size();
	};
	Matrix dot(const Matrix& m) const;
	std::vector<double> dot(const std::vector<double>& v) const;
	void print() const;

	static Matrix eye(size_t n, double a = 1);
};

std::vector<double> solveGaussFullPivot(
	const Matrix& A,
	const std::vector<double>& b,
	double eps = 1e-12
);

class FEM {
	size_t mx, ny;
	Point left_down;
	Point right_up;
	std::vector<Point> points;
	std::vector<Triangle> triangles;
	double triangle_area;
	std::vector<Point> f;
	std::vector<Point> u;
	CoordinateSystem coordinate_system;

	std::vector<double> F;
	Matrix A;

public:
	FEM(const Point& a, const Point& b, size_t n, size_t m,
		CoordinateSystem coordinate_system = CoordinateSystem::Axisymmetric);

	Point& get_point(size_t n);
	const Point& get_point(size_t n) const;
	Triangle& get_triangle(size_t n);
	const Triangle& get_triangle(size_t n) const;

	size_t psize() const { return points.size(); };
	size_t tsize() const { return triangles.size(); };
	size_t xsize() const { return mx; };
	size_t ysize() const { return ny; };
	CoordinateSystem get_coordinate_system() const { return coordinate_system; };

	void print_points() const;
	void print_triangles() const;

	Point& operator[](size_t n) { return get_point(n); };
	const Point& operator[](size_t n) const { return get_point(n); };

	Triangle& operator()(size_t n) { return get_triangle(n); };
	const Triangle& operator()(size_t n) const { return get_triangle(n); };

	void set_boundaries(char side, const vec_function& g);

	void construct_AF(double E, double nu, vec_function body_force);
	void apply_boundaries();
	std::vector<Point> solve();

	void clear_AFu();

	void bc2_side(lambda_func j, size_t start, size_t finish, double len,
		int side, const std::vector<vec_function>& g,
		std::vector<double>& p_vec);

	void calculate_bc2(const std::vector<size_t>& pos,
		const std::vector<vec_function>& g, std::vector<double>& p_vec);

	std::pair<Matrix, std::vector<double>> get_AF();
	void set_AF(const Matrix& A_new, const std::vector<double>& F_new);
};

Matrix operator*(double a, Matrix m);

std::tuple<Matrix, Matrix> LU_decomposition(const Matrix& m);

std::vector<double> solveLU(const Matrix& L, const Matrix& U,
	const std::vector<double>& b);

Point zero(const Point& p);

class FSEM {
	double E;
	double nu;
	Point a;
	Point b;
	size_t n_side_x;
	size_t n_side_y;
	std::vector<Point> nodes;
	size_t coef_x;
	size_t coef_y;

	std::vector<std::vector<Point>> basis;
	std::vector<Point> basis_coefficients;

	Matrix K;
	std::vector<double> f;
	CoordinateSystem coordinate_system;

public:
	FEM fem;

	FSEM(double E, double nu, const Point& a, const Point& b,
		size_t n_x, size_t n_y, int coef_val_x = 1, int coef_val_y = 1,
		CoordinateSystem coordinate_system = CoordinateSystem::Axisymmetric);

	Point& get_node(size_t n);
	const Point& get_node(size_t n) const;

	size_t nsize() const { return nodes.size(); };
	CoordinateSystem get_coordinate_system() const { return coordinate_system; };

	Point& operator[](size_t n) { return get_node(n); };
	const Point& operator[](size_t n) const { return get_node(n); };

	void print_nodes() const;

	void construct_basis();

	Matrix matrix_form_basis();

	const Matrix& get_K() const { return K; }
	const std::vector<double>& get_f() const { return f; }
	const std::vector<std::vector<Point>>& get_basis() const { return basis; }

	void construct_f_bc2(const std::vector<size_t>& pos,
		const std::vector<vec_function>& g);

	std::vector<Point> find_answer();
	std::vector<Point> find_answer(const std::vector<double>& coefs, size_t start = 0);

	void set_bc1(char side, const vec_function& g);

	void calculate_coef_Matrix_bc2(const int finish, const int i,
		bool cur_pos, bool prev_pos, int& add_B, int& add_C, const int add_basis,
		Matrix& B, Matrix& C);

	void save_bc1(std::vector<double>& coefs_Dirichle, int& dir_id, const int i);

	void set_bc2(const std::vector<size_t>& pos,
		const std::vector<vec_function>& g);

	std::vector<size_t> get_side_nodes(char side) const;
	std::vector<size_t> get_side_fem_nodes(char side) const;

	Point coefficient(int i, Point coef_val) {
		return std::isnan(basis_coefficients[i].x) ? coef_val : basis_coefficients[i];
	}

	std::vector<std::pair<size_t, double>> get_known_dofs() const;
};

double mortar_shape_func(size_t i, const std::vector<double>& s, double cur);

std::vector<double> solve_mortar_contact(
	FSEM& bottom_body,
	FSEM& top_body,
	const std::vector<double>& rhs_bottom,
	const std::vector<double>& rhs_top,
	const ContactOptions& options);

std::vector<double> solve_mortar_contact(
	FSEM& bottom_body,
	FSEM& top_body,
	const std::vector<double>& rhs_bottom,
	const std::vector<double>& rhs_top,
	size_t lambda_node_count = 0);

std::vector<double> solveWithLU(const Matrix& A,
	const std::vector<double>& b,
	double eps = 1e-15);
