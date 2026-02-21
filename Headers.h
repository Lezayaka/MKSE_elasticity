#pragma once
#include <vector>
#include <iostream>
#include <functional>

struct Point {
	double x, y;

	friend std::ostream& operator<<(std::ostream& output, const Point& p);
};

Point operator+(const Point& a, const Point& b);
Point operator-(const Point& a, const Point& b);
Point operator*(const double a, const Point& p);
Point operator/(const Point& p, const double a);

std::vector<double> operator-(const std::vector<double>& a,
	const std::vector<double>& b);

// Вукции вида f(x, y) = z
using function = std::function<double(const Point&)>;
// Вукции вида f(x, y) = {u, v}
using vec_function = std::function<Point(const Point&)>;

using lambda_func = std::function<int(int)>;

struct Triangle {
	size_t a, b, c;
};

struct Rectangle {
	size_t a, b, c, d;
};

// Класс для работы с матрицами
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

// Класс хранящий сетку для МКЭ
class FEM {
	size_t mx, ny;
	Point left_down;
	Point right_up;
	std::vector<Point> points;
	std::vector<Triangle> triangles;
	double triangle_area;
	std::vector<Point> f;
	std::vector<Point> u;
	
	std::vector<double> F;
	Matrix A;

public:
	// значение решения в узле
	FEM(const Point& a, const Point& b, size_t n, size_t m);

	Point& get_point(size_t n);
	const Point& get_point(size_t n) const;
	Triangle& get_triangle(size_t n);
	const Triangle& get_triangle(size_t n) const;

	size_t psize() const { return points.size(); };
	size_t tsize() const { return triangles.size(); };

	void print_points() const;
	void print_triangles() const;

	// Получение точки сетки
	Point& operator[](size_t n) { return get_point(n); };
	const Point& operator[](size_t n) const { return get_point(n); };

	// Получение элемента сетки
	Triangle& operator()(size_t n) { return get_triangle(n); };
	const Triangle& operator()(size_t n) const {
		return get_triangle(n);
	};

	// +-N-+ ♡♡♡♡♡♡♡♡♡♡♡♡
	// W---E ♡♡♡♡♡♡♡♡♡♡♡♡
	// +-S-+ ♡♡♡♡♡♡♡♡♡♡♡♡
	void set_boundaries(char side, const vec_function& g);

	void construct_AF(double E, double nu, vec_function f);
	void apply_boundaries();
	std::vector<Point> solve();

	void clear_AFu();

	void bc2_side(lambda_func j, int start, int finish, double len, 
		int side, bool prev_side, bool next_side,
		const std::vector<vec_function>& g, std::vector<double>& p_vec);

	void calculate_bc2(const std::vector<size_t>& pos,
		const std::vector<vec_function>& g, std::vector<double>& p_vec);

	std::pair<Matrix, std::vector<double>> get_AF();
	void set_AF(const Matrix& A_new, const std::vector<double>& F_new);
};

Matrix operator*(double a, Matrix m);

std::tuple<Matrix, Matrix> LU_decomposition(const Matrix& m);

std::vector<double> solveLU(const Matrix& L, const Matrix& U, const std::vector<double>& b);

Point zero(const Point& p);

// Класс для метода конечных суперэлементов
class FSEM {
	double E;
	double nu;
	Point a;
	Point b;
	size_t n_side_x; // количество отрезков на горизонтальной границе 
	size_t n_side_y; // количество отрезков на вертикальной границе
	// координаты узлов, начиная с левого нижнего, идут по часовой стрелке
	std::vector<Point> nodes; 
	//int coef_x; // количество узлов сетки для МКЭ по x = n_side_x * coef
	//int coef_y; // количество узлов сетки для МКЭ по y = n_side_y * coef

	// храним сеточные значения базисных функций
	std::vector<std::vector<Point>> basis;

	// коэффициенты
	std::vector<Point> basis_coefficients;

	Matrix K; // матрица жесткости
	std::vector<double> f;

public:
	FEM fem;

	FSEM(double E, double nu, const Point& a, const Point& b, 
		size_t n_x, size_t n_y, int coef_val_x = 1, int coef_val_y = 1);

	Point& get_node(size_t n);
	const Point& get_node(size_t n) const;

	size_t nsize() const { return nodes.size(); };

	// получение узла сетки
	Point& operator[](size_t n) { return get_node(n); };
	const Point& operator[](size_t n) const { return get_node(n); };

	void print_nodes() const;

	void construct_basis();

	Matrix matrix_form_basis();

	// интегралы от ГУ 2 рода * функции формы мкэ
	void construct_f_bc2(const std::vector<size_t>& pos,
		const std::vector<vec_function>& g);

	//void find_coefficients(const vec_function& g);
	std::vector<Point> find_answer();

	// +-N-+ ♡♡♡♡♡♡♡♡♡♡♡♡
	// W---E ♡♡♡♡♡♡♡♡♡♡♡♡
	// +-S-+ ♡♡♡♡♡♡♡♡♡♡♡♡
	void set_bc1(char side, const vec_function& g);

	void calculate_coef_Matrix_bc2(const int finish, const int i,
		bool cur_pos, bool prev_pos, int& add_B, int& add_C, const int add_basis, Matrix& B,
		Matrix& C);

	void save_bc1(std::vector<double>& coefs_Dirichle, int& dir_id, const int i);

	// pos = {W, N, E, S} 
	void set_bc2(const std::vector<size_t>& pos, 
		const std::vector<vec_function>& g);
};