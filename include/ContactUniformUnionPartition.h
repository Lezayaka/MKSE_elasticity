#pragma once

// Contact discretization on a uniform refinement of the combined traces.

#include <algorithm>
#include <cmath>

#include "Elasticity.h"
#include "ContactUniformLambdaPartition.h"

namespace contact_uniform_union_partition {

constexpr double kContactEps = 1e-12;

inline bool almost_equal(double lhs, double rhs) {
	return std::fabs(lhs - rhs) < kContactEps;
}

inline std::vector<MortarElement> build_mortar_elements(
	const std::vector<double>& bottom_x,
	const std::vector<double>& top_x,
	const std::vector<double>& lambda_nodes,
	double contact_left,
	double contact_right) {

	std::vector<double> partition;
	partition.reserve(bottom_x.size() + top_x.size() + lambda_nodes.size() + 2);
	partition.push_back(contact_left);
	partition.push_back(contact_right);

	auto append_contact_nodes = [&](const std::vector<double>& nodes) {
		for (double x : nodes)
			if (x >= contact_left - kContactEps && x <= contact_right + kContactEps)
				partition.push_back(x);
	};

	append_contact_nodes(bottom_x);
	append_contact_nodes(top_x);
	append_contact_nodes(lambda_nodes);

	std::sort(partition.begin(), partition.end());
	partition.erase(std::unique(partition.begin(), partition.end(), almost_equal), partition.end());

	std::vector<MortarElement> mortar_elements;
	mortar_elements.reserve(partition.size() - 1);
	for (size_t seg = 0; seg + 1 < partition.size(); ++seg)
		mortar_elements.push_back({ partition[seg], partition[seg + 1] });

	return mortar_elements;
}

inline ContactDiscretization build(
	const std::vector<double>& bottom_x,
	const std::vector<double>& top_x,
	double contact_left,
	double contact_right,
	size_t lambda_node_count) {

	ContactDiscretization discretization =
		contact_uniform_lambda_partition::build(
			bottom_x, top_x, contact_left, contact_right, lambda_node_count);

	discretization.mortar_elements = build_mortar_elements(
		bottom_x,
		top_x,
		discretization.lambda_nodes,
		contact_left,
		contact_right);

	return discretization;
}

} // namespace contact_uniform_union_partition
