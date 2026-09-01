#pragma once

// Contact discretization on an independent uniform multiplier grid.

#include <algorithm>
#include <cmath>

#include "Elasticity.h"

namespace contact_uniform_lambda_partition {

constexpr double kContactEps = 1e-12;

inline size_t count_contact_nodes(
	const std::vector<double>& x_nodes,
	double contact_left,
	double contact_right) {

	size_t count = 0;
	for (double x : x_nodes)
		if (x >= contact_left - kContactEps && x <= contact_right + kContactEps)
			++count;

	return count;
}

inline size_t default_lambda_node_count(
	const std::vector<double>& bottom_x,
	const std::vector<double>& top_x,
	double contact_left,
	double contact_right) {

	return std::max<size_t>(
		2,
		std::max(
			count_contact_nodes(bottom_x, contact_left, contact_right),
			count_contact_nodes(top_x, contact_left, contact_right)));
}

inline std::vector<double> build_uniform_lambda_nodes(
	double contact_left,
	double contact_right,
	size_t lambda_node_count) {

	std::vector<double> lambda_nodes(lambda_node_count);
	const double step = (contact_right - contact_left) / (lambda_node_count - 1);

	for (size_t i = 0; i < lambda_node_count; ++i)
		lambda_nodes[i] = contact_left + i * step;

	lambda_nodes.front() = contact_left;
	lambda_nodes.back() = contact_right;
	return lambda_nodes;
}

inline ContactDiscretization build(
	const std::vector<double>& bottom_x,
	const std::vector<double>& top_x,
	double contact_left,
	double contact_right,
	size_t lambda_node_count) {

	if (lambda_node_count == 0)
		lambda_node_count = default_lambda_node_count(
			bottom_x, top_x, contact_left, contact_right);

	ContactDiscretization discretization;
	discretization.lambda_nodes = build_uniform_lambda_nodes(
		contact_left, contact_right, std::max<size_t>(2, lambda_node_count));

	discretization.mortar_elements.reserve(discretization.lambda_nodes.size() - 1);
	for (size_t seg = 0; seg + 1 < discretization.lambda_nodes.size(); ++seg)
		discretization.mortar_elements.push_back({
			discretization.lambda_nodes[seg],
			discretization.lambda_nodes[seg + 1]
		});

	return discretization;
}

} // namespace contact_uniform_lambda_partition
