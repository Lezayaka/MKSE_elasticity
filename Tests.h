#pragma once

#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <string>

#include "Headers.h"

struct MeshSize {
	size_t bottom_x;
	size_t bottom_y;
	size_t top_x;
	size_t top_y;
	size_t lambda_nodes;
};

struct TestCase {
	std::function<vec_function(double, double)> exact_solution;
	double E;
	double nu;
	Point bottom_a;
	Point bottom_b;
	Point top_a;
	Point top_b;
	MeshSize mesh;
	std::array<char, 3> bottom_dirichlet_sides;
	std::array<char, 3> top_dirichlet_sides;
};

const std::map<std::string, TestCase> TESTS = {
	{
		"axisymmetric_inverse_r",
		TestCase{
			[](double, double) {
				return [](const Point& p) {
					return Point{ 5.0 / p.x, 0.0 };
				};
			},
			21e10,
			0.3,
			Point{ 1.0, 0.0 },
			Point{ 3.0, 0.5 },
			Point{ 1.0, 0.5 },
			Point{ 3.0, 3.0 },
			MeshSize{ 18, 18, 18, 18, 18 },
			std::array<char, 3>{ 'W', 'E', 'S' },
			std::array<char, 3>{ 'W', 'N', 'E' }
		}
	},
	{
		"axisymmetric_linear",
		TestCase{
			[](double, double) {
				return [](const Point& p) {
					return Point{ 2.0 * p.x + 5.0 / p.x, 34.0 * p.y };
				};
			},
			21e10,
			0.3,
			Point{ 1.0, 0.0 },
			Point{ 3.0, 0.5 },
			Point{ 1.0, 0.5 },
			Point{ 3.0, 3.0 },
			MeshSize{ 10, 10, 10, 10, 10 },
			std::array<char, 3>{ 'W', 'E', 'S' },
			std::array<char, 3>{ 'W', 'N', 'E' }
		}
	},
	{
		"cartesian_exp",
		TestCase{
			[](double, double) {
				return [](const Point& p) {
					return Point{
						std::exp(p.x) * std::cos(p.y - 0.5),
						-std::exp(p.x) * std::sin(p.y - 0.5)
					};
				};
			},
			21e10,
			0.3,
			Point{ 0.0, 0.0 },
			Point{ 1.0, 0.5 },
			Point{ 0.0, 0.5 },
			Point{ 1.0, 1.0 },
			MeshSize{ 6, 6, 10, 10, 10 },
			std::array<char, 3>{ 'W', 'E', 'S' },
			std::array<char, 3>{ 'W', 'N', 'E' }
		}
	},
	{
		"cartesian_linear",
		TestCase{
			[](double, double) {
				return [](const Point& p) {
					return Point{ -21.0 * p.x, 13.0 * p.y };
				};
			},
			21e10,
			0.3,
			Point{ 0.0, 0.0 },
			Point{ 3.0, 1.0 },
			Point{ 0.0, 1.0 },
			Point{ 2.0, 4.0 },
			MeshSize{ 5, 5, 10, 10, 10 },
			std::array<char, 3>{ 'W', 'E', 'S' },
			std::array<char, 3>{ 'W', 'N', 'E' }
		}
	}
};
