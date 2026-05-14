#pragma once

#include "Headers.h"

namespace fem_cartesian {

inline void assemble(
	Matrix& A,
	std::vector<double>& F,
	const std::vector<Point>& points,
	const std::vector<Triangle>& triangles,
	double triangle_area,
	double E,
	double nu,
	const vec_function& f) {

	const double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
	const double mu = E / (2.0 * (1.0 + nu));
	const Matrix C({
		{ lambda + 2.0 * mu, lambda, 0.0 },
		{ lambda, lambda + 2.0 * mu, 0.0 },
		{ 0.0, 0.0, mu }
	});

	for (const Triangle& T : triangles) {
		const Point& p1 = points[T.a];
		const Point& p2 = points[T.b];
		const Point& p3 = points[T.c];
		const Matrix grad_phi({
			{ (p2.y - p3.y) / (2.0 * triangle_area), (p3.x - p2.x) / (2.0 * triangle_area) },
			{ (p3.y - p1.y) / (2.0 * triangle_area), (p1.x - p3.x) / (2.0 * triangle_area) },
			{ (p1.y - p2.y) / (2.0 * triangle_area), (p2.x - p1.x) / (2.0 * triangle_area) }
		});

		Matrix R(3, 6ull);
		for (size_t p = 0; p < 3; ++p) {
			R[0][2 * p] = grad_phi[p][0];
			R[1][2 * p + 1] = grad_phi[p][1];
			R[2][2 * p] = grad_phi[p][1];
			R[2][2 * p + 1] = grad_phi[p][0];
		}

		Matrix Ae = triangle_area * R.T().dot(C).dot(R);

		std::vector<double> Fe(6, 0.0);
		const Point center = (p1 + p2 + p3) / 3.0;
		const Point body_force = f(center);
		for (size_t p = 0; p < 3; ++p) {
			Fe[2 * p] = body_force.x * triangle_area / 3.0;
			Fe[2 * p + 1] = body_force.y * triangle_area / 3.0;
		}

		auto dof = [T](size_t i) -> size_t {
			if (i == 0) return 2 * T.a;
			if (i == 1) return 2 * T.a + 1;
			if (i == 2) return 2 * T.b;
			if (i == 3) return 2 * T.b + 1;
			if (i == 4) return 2 * T.c;
			return 2 * T.c + 1;
		};

		for (size_t i = 0; i < 6; ++i) {
			F[dof(i)] += Fe[i];
			for (size_t j = 0; j < 6; ++j)
				A[dof(i)][dof(j)] += Ae[i][j];
		}
	}
}

inline double boundary_segment_weight(const Point&, double len) {
	return 0.5 * len;
}

} // namespace fem_cartesian
