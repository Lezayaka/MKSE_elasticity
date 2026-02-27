#include "Headers.h"

Matrix::Matrix(size_t n, double a) {
	matrix = std::vector<std::vector<double>>(n, std::vector<double>(n, a));
}

Matrix::Matrix(size_t n, size_t m, double a) {
	matrix = std::vector<std::vector<double>>(n, std::vector<double>(m, a));
}

Matrix Matrix::operator*(double a) {
	Matrix m = *this;

	for (size_t i = 0; i != matrix.size(); ++i) {
		for (size_t j = 0; j != matrix[i].size(); ++j) {
			m[i][j] *= a;
		}
	}

	return m;
}

Matrix operator*(double a, Matrix m) {
	return m * a;
}

Matrix& Matrix::operator*=(double a) {
	for (size_t i = 0; i != matrix.size(); ++i) {
		for (size_t j = 0; j != matrix[i].size(); ++j) {
			matrix[i][j] *= a;
		}
	}

	return *this;
}

std::vector<double> operator-(const std::vector<double>& a, 
	const std::vector<double>& b) {
	if (a.size() != b.size()) 
		std::cout << "Wrong vector size for subtraction\n";

	std::vector<double> res(a.size());

	for (int i = 0; i < a.size(); ++i)
		res[i] = a[i] - b[i];

	return res;
}

Matrix Matrix::T() const {
	Matrix m(size(1), size());

	for (size_t i = 0; i != size(1); ++i) {
		for (size_t j = 0; j != size(); ++j) {
			m[i][j] = matrix[j][i];
		}
	}

	return m;
}

Matrix Matrix::dot(const Matrix& m) const {
	Matrix res(size(), m.size(1));

	if (m.size() != size(1))
		throw std::runtime_error("Wrong matrix size while multiplication!");

	for (size_t i = 0; i != res.size(); ++i) {
		for (size_t j = 0; j != res.size(1); ++j) {

			for (size_t k = 0; k != m.size(); ++k) {
				res[i][j] += matrix[i][k] * m[k][j];
			}
		}
	}

	return res;
}

std::vector<double> Matrix::dot(const std::vector<double>&v) const {
	if (v.size() != size(1))
		throw std::runtime_error("Wrong matrix size while multiplication!");

	std::vector<double> res(size());

	for (int i = 0; i < size(); ++i)
		for (int j = 0; j < size(1); ++j)
			res[i] += matrix[i][j] * v[j];

	return res;
}

void Matrix::print() const {
	for (size_t i = 0; i != size(); ++i) {
		std::cout << "{";
		for (size_t j = 0; j != size(1); ++j) {
			std::cout << matrix[i][j] << ", ";
		}
		std::cout << "},\n";
	}
}

Matrix Matrix::eye(size_t n, double a) {
	Matrix m(n, n);
	for (size_t i = 0; i != m.size(); ++i) {
		m[i][i] = a;
	}

	return m;
}

std::tuple<Matrix, Matrix> LU_decomposition(const Matrix& m) {
	Matrix U(m.size(), m.size(1));
	Matrix L = Matrix::eye(m.size());

	for (size_t i = 0; i != m.size(); ++i) {
		for (size_t j = 0; j != m.size(); ++j) {
			if (i <= j) {
				U[i][j] = m[i][j];
				for (size_t k = 0; k != i; ++k) {
					U[i][j] -= L[i][k] * U[k][j];
				}
			}

			else {
				L[i][j] = m[i][j];
				for (size_t k = 0; k != j; ++k) {
					L[i][j] -= L[i][k] * U[k][j];
				}
				L[i][j] /= U[j][j];
			}
		}
	}

	return { L, U };
}

std::vector<double> solveLU(const Matrix& L, const Matrix& U,
	const std::vector<double>& b) {
	std::vector<double> y(L.size());
	for (size_t i = 0; i != L.size(); ++i) {
		y[i] = b[i];
		for (size_t j = 0; j != i; ++j) {
			y[i] -= L[i][j] * y[j];
		}
	}

	std::vector<double> x(L.size());
	for (size_t i = L.size() - 1; i + 1 != 0; --i) {
		x[i] = y[i];
		for (size_t j = i + 1; j != L.size(); ++j) {
			x[i] -= U[i][j] * x[j];
		}
		x[i] /= U[i][i];
	}

	return x;
}

std::vector<double> solveGaussFullPivot(const Matrix& A,
	const std::vector<double>& b,
	double eps) {
	
	const size_t n = A.size();
	Matrix M(A);
	std::vector<double> rhs = b;

	// перестановки столбцов
	std::vector<size_t> col_perm(n);

	for (size_t i = 0; i != n; ++i)
		col_perm[i] = i;

	for (size_t k = 0; k != n; ++k) {
		size_t pivot_row = k;
		size_t pivot_col = k;
		double pivot_abs = 0;

		for (size_t i = k; i != n; ++i) {
			for (size_t j = k; j != n; ++j) {
				double cur = abs(M[i][j]);
				if (cur > pivot_abs) {
					pivot_abs = cur;
					pivot_row = i;
					pivot_col = j;
				}
			}
		}

		if (pivot_row != k) {
			std::swap(M[pivot_row], M[k]);
			std::swap(rhs[pivot_row], rhs[k]);
		}

		if (pivot_col != k) {
			for (size_t i = 0; i != n; ++i)
				std::swap(M[i][pivot_col], M[i][k]);
			std::swap(col_perm[pivot_col], col_perm[k]);
		}

		for (size_t i = k + 1; i != n; ++i) {
			double factor = M[i][k] / M[k][k];
			M[i][k] = 0;
			for (size_t j = k + 1; j != n; ++j)
				M[i][j] -= factor * M[k][j];
			rhs[i] -= factor * rhs[k];
		}
	}

	std::vector<double> y(n);
	for (size_t i = n; i-- > 0;) {
		double sum = rhs[i];
		for (size_t j = i + 1; j != n; ++j)
			sum -= M[i][j] * y[j];

		y[i] = sum / M[i][i];
	}

	std::vector<double> x(n);
	for (size_t i = 0; i != n; ++i)
		x[col_perm[i]] = y[i];

	return x;
}