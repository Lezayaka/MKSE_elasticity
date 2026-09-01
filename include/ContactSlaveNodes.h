#pragma once

// Contact discretization based on the passive body's trace nodes.

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "Elasticity.h"

namespace contact_slave_nodes {

constexpr double kContactEps = 1e-12;

inline bool almost_equal(double lhs, double rhs) {
	return std::fabs(lhs - rhs) < kContactEps;
}

inline std::vector<double> collect_overlap_nodes(
	const std::vector<double>& bottom_x,
	const std::vector<double>& top_x,
	double contact_left,
	double contact_right) {

	std::vector<double> nodes;
	nodes.reserve(bottom_x.size() + top_x.size() + 2);
	nodes.push_back(contact_left);
	nodes.push_back(contact_right);

	auto append_inside = [&](const std::vector<double>& x_nodes) {
		for (double x : x_nodes)
			if (x >= contact_left - kContactEps && x <= contact_right + kContactEps)
				nodes.push_back(x);
	};

	append_inside(bottom_x);
	append_inside(top_x);

	std::sort(nodes.begin(), nodes.end());
	nodes.erase(std::unique(nodes.begin(), nodes.end(), almost_equal), nodes.end());
	return nodes;
}

inline std::vector<double> collect_active_lambda_nodes(
	const std::vector<double>& passive_x,
	double contact_left,
	double contact_right) {

	std::vector<double> active_nodes;
	active_nodes.reserve(passive_x.size());

	for (size_t i = 0; i < passive_x.size(); ++i) {
		const double support_left = (i == 0) ? passive_x[i] : passive_x[i - 1];
		const double support_right = (i + 1 == passive_x.size()) ? passive_x[i] : passive_x[i + 1];

		if (support_right > contact_left + kContactEps &&
			support_left < contact_right - kContactEps) {
			active_nodes.push_back(passive_x[i]);
		}
	}

	if (active_nodes.size() < 2)
		throw std::runtime_error("Not enough passive contact nodes for Lagrange multipliers.");

	return active_nodes;
}

inline ContactDiscretization build(
	const std::vector<double>& bottom_x,
	const std::vector<double>& top_x,
	double contact_left,
	double contact_right,
	ContactSlaveBody passive_body) {

	const std::vector<double>& passive_x =
		(passive_body == ContactSlaveBody::Bottom) ? bottom_x : top_x;

	ContactDiscretization discretization;
	discretization.lambda_nodes = collect_active_lambda_nodes(
		passive_x, contact_left, contact_right);

	const std::vector<double> integration_nodes = collect_overlap_nodes(
		bottom_x, top_x, contact_left, contact_right);
	discretization.mortar_elements.reserve(integration_nodes.size() - 1);
	for (size_t seg = 0; seg + 1 < integration_nodes.size(); ++seg)
		discretization.mortar_elements.push_back({ integration_nodes[seg], integration_nodes[seg + 1] });

	return discretization;
}

} // namespace contact_slave_nodes
