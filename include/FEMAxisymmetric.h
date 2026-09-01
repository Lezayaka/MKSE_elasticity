#pragma once

// Axisymmetric finite-element integration kernels.

#include <array>
#include <cmath>

#include "Elasticity.h"

namespace fem_axisymmetric {

constexpr double kTwoPi = 6.28318530717958647692;

struct TriangleQuadraturePoint {
	double weight;
	std::array<double, 3> phi;
};

inline const std::array<TriangleQuadraturePoint, 3> kTriangleQuadrature = { {
	{ 1.0 / 3.0, { 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0 } },
	{ 1.0 / 3.0, { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 } },
	{ 1.0 / 3.0, { 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0 } }
} };

inline double signed_double_area(const Point& p1, const Point& p2, const Point& p3) {
	return (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y);
}

inline void assemble(
	Matrix& A,
	std::vector<double>& F,
	const std::vector<Point>& points,
	const std::vector<Triangle>& triangles,
	double,
	double E,
	double nu,
	const vec_function& f) {

	const double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
	const double mu = E / (2.0 * (1.0 + nu));
	const Matrix C({
		{ lambda + 2.0 * mu, lambda, lambda, 0.0 },
		{ lambda, lambda + 2.0 * mu, lambda, 0.0 },
		{ lambda, lambda, lambda + 2.0 * mu, 0.0 },
		{ 0.0, 0.0, 0.0, mu }
	});

	for (const Triangle& T : triangles) {
		const Point& p1 = points[T.a];
		const Point& p2 = points[T.b];
		const Point& p3 = points[T.c];
		const double two_area = signed_double_area(p1, p2, p3);
		const double area = 0.5 * std::fabs(two_area);

		const Matrix grad_phi({
			{ (p2.y - p3.y) / two_area, (p3.x - p2.x) / two_area },
			{ (p3.y - p1.y) / two_area, (p1.x - p3.x) / two_area },
			{ (p1.y - p2.y) / two_area, (p2.x - p1.x) / two_area }
		});

		Matrix Ae(6ull);
		std::vector<double> Fe(6, 0.0);
		for (const auto& qp : kTriangleQuadrature) {
			const Point q = qp.phi[0] * p1 + qp.phi[1] * p2 + qp.phi[2] * p3;
			const double r_q = q.x;

			Matrix B(4ull, 6ull);
			for (size_t p = 0; p < 3; ++p) {
				const double dphi_dr = grad_phi[p][0];
				const double dphi_dz = grad_phi[p][1];
				const double phi = qp.phi[p];

				B[0][2 * p] = dphi_dr;
				B[1][2 * p + 1] = dphi_dz;
				B[2][2 * p] = phi / r_q;
				B[3][2 * p] = dphi_dz;
				B[3][2 * p + 1] = dphi_dr;
			}

			const Matrix stiffness_q = B.T().dot(C).dot(B);
			const double weight = kTwoPi * area * qp.weight * r_q;
			for (size_t i = 0; i < 6; ++i)
				for (size_t j = 0; j < 6; ++j)
					Ae[i][j] += weight * stiffness_q[i][j];

			const Point body_force = f(q);
			for (size_t p = 0; p < 3; ++p) {
				const double shape_value = qp.phi[p];
				Fe[2 * p] += weight * shape_value * body_force.x;
				Fe[2 * p + 1] += weight * shape_value * body_force.y;
			}
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

inline double boundary_segment_weight(const Point& midpoint, double len) {
	return 0.5 * kTwoPi * midpoint.x * len;
}

} // namespace fem_axisymmetric
