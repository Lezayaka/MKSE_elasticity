#include "Headers.h"
#include "Tests.h"
#include "ContactSlaveNodes.h"
#include "ContactUniformLambdaPartition.h"
#include "ContactUniformUnionPartition.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace {

constexpr double kTraceEps = 1e-12;

struct TestRunConfig {
	std::string name;
	std::optional<size_t> bottom_x;
	std::optional<size_t> bottom_y;
	std::optional<size_t> top_x;
	std::optional<size_t> top_y;
	std::optional<size_t> lambda_nodes;
};

struct InputConfig {
	CoordinateSystem coordinate_system = CoordinateSystem::Axisymmetric;
	ContactMethod contact_method = ContactMethod::UniformUnionPartition;
	ContactSlaveBody slave_body = ContactSlaveBody::Bottom;
	std::vector<TestRunConfig> tests;
};

std::string trim(const std::string& text) {
	std::string value = text;
	if (value.size() >= 3 &&
		static_cast<unsigned char>(value[0]) == 0xEF &&
		static_cast<unsigned char>(value[1]) == 0xBB &&
		static_cast<unsigned char>(value[2]) == 0xBF) {
		value.erase(0, 3);
	}

	const auto first = value.find_first_not_of(" \t\r\n");
	if (first == std::string::npos)
		return {};
	const auto last = value.find_last_not_of(" \t\r\n");
	return value.substr(first, last - first + 1);
}

std::string normalize(std::string value) {
	value = trim(value);
	if (value.size() >= 2 &&
		((value.front() == '"' && value.back() == '"') ||
			(value.front() == '\'' && value.back() == '\''))) {
		value = value.substr(1, value.size() - 2);
	}
	std::replace(value.begin(), value.end(), '-', '_');
	for (char& ch : value)
		ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
	return value;
}

std::vector<std::string> split_list(std::string value) {
	for (char& ch : value)
		if (ch == ',' || ch == ';')
			ch = ' ';

	std::vector<std::string> result;
	std::istringstream in(value);
	std::string token;
	while (in >> token)
		result.push_back(token);
	return result;
}

bool starts_with(const std::string& value, const std::string& prefix) {
	return value.size() >= prefix.size() &&
		std::equal(prefix.begin(), prefix.end(), value.begin());
}

size_t parse_size(const std::string& value, const std::string& key) {
	size_t parsed = 0;
	size_t used = 0;
	parsed = std::stoull(value, &used);
	if (used != value.size())
		throw std::runtime_error("Invalid integer value for '" + key + "': " + value);
	return parsed;
}

CoordinateSystem parse_coordinate_system(const std::string& value) {
	const std::string v = normalize(value);
	if (v == "a" || v == "axi" || v == "axisymmetric" || v == "axisym" ||
		v == "r_z" || v == "осесимметричная" || v == "осесимметричной" ||
		v == "осисиметричная" || v == "осисиметричной") {
		return CoordinateSystem::Axisymmetric;
	}
	if (v == "c" || v == "cartesian" || v == "cart" || v == "x_y" ||
		v == "dsk" || v == "дск" || trim(value) == "ДСК" ||
		v == "декартовая" || v == "декартовой") {
		return CoordinateSystem::Cartesian;
	}
	throw std::runtime_error("Unknown coordinate system: " + value);
}

ContactMethod parse_contact_method(const std::string& value) {
	const std::string v = normalize(value);
	if (v == "s" || v == "slave" || v == "slave_nodes" ||
		v == "passive" || v == "узлы" || v == "пассивное") {
		return ContactMethod::SlaveNodes;
	}
	if (v == "d" || v == "lambda" || v == "uniform_lambda" ||
		v == "uniform" || v == "dissertation" || v == "диссертация") {
		return ContactMethod::UniformLambdaPartition;
	}
	if (v == "u" || v == "all" || v == "union" || v == "uniform_union" ||
		v == "po_vsem" || v == "по_всем") {
		return ContactMethod::UniformUnionPartition;
	}
	throw std::runtime_error("Unknown contact method: " + value);
}

ContactSlaveBody parse_slave_body(const std::string& value) {
	const std::string v = normalize(value);
	if (v == "bottom" || v == "b" || v == "lower" || v == "низ" || v == "нижнее")
		return ContactSlaveBody::Bottom;
	if (v == "top" || v == "t" || v == "upper" || v == "верх" || v == "верхнее")
		return ContactSlaveBody::Top;
	throw std::runtime_error("Unknown passive body: " + value);
}

TestRunConfig& append_test(InputConfig& config, const std::string& name) {
	config.tests.push_back(TestRunConfig{ name });
	return config.tests.back();
}

void add_test_name(InputConfig& config, const std::string& name) {
	if (normalize(name) == "all") {
		for (const auto& [test_name, test] : TESTS) {
			(void)test;
			append_test(config, test_name);
		}
	}
	else {
		append_test(config, name);
	}
}

void apply_override(TestRunConfig& run, const std::string& key, const std::string& value) {
	const std::string k = normalize(key);
	const std::string v = normalize(value);
	if (k == "bottom_x" || k == "bottom_nx" || k == "n_bottom_x" || k == "bx")
		run.bottom_x = parse_size(v, key);
	else if (k == "bottom_y" || k == "bottom_ny" || k == "n_bottom_y" || k == "by")
		run.bottom_y = parse_size(v, key);
	else if (k == "top_x" || k == "top_nx" || k == "n_top_x" || k == "tx")
		run.top_x = parse_size(v, key);
	else if (k == "top_y" || k == "top_ny" || k == "n_top_y" || k == "ty")
		run.top_y = parse_size(v, key);
	else if (k == "lambda" || k == "lambda_nodes" || k == "n_lambda" || k == "multipliers")
		run.lambda_nodes = parse_size(v, key);
	else
		throw std::runtime_error("Unknown per-test parameter: " + key);
}

void apply_key_value(InputConfig& config, TestRunConfig* current_test,
	const std::string& key, const std::string& value) {

	const std::string k = normalize(key);
	if (k == "coordinate" || k == "coordinates" || k == "coord" ||
		k == "system" || k == "система") {
		config.coordinate_system = parse_coordinate_system(value);
	}
	else if (k == "contact" || k == "method" || k == "contact_method" || k == "контакт") {
		config.contact_method = parse_contact_method(value);
	}
	else if (k == "passive" || k == "slave" || k == "slave_body" || k == "passive_body") {
		config.slave_body = parse_slave_body(value);
	}
	else if (k == "tests" || k == "test" || k == "тесты") {
		for (const std::string& name : split_list(value))
			add_test_name(config, name);
	}
	else if (current_test != nullptr) {
		apply_override(*current_test, key, value);
	}
	else {
		throw std::runtime_error("Unknown input key outside a test block: " + key);
	}
}

void parse_assignment_line(InputConfig& config, TestRunConfig* current_test, const std::string& line) {
	const size_t eq = line.find('=');
	const size_t colon = line.find(':');
	const size_t pos = std::min(
		eq == std::string::npos ? line.size() : eq,
		colon == std::string::npos ? line.size() : colon);

	if (pos != line.size()) {
		apply_key_value(config, current_test, line.substr(0, pos), line.substr(pos + 1));
		return;
	}

	std::istringstream in(line);
	std::string key;
	std::string value;
	in >> key >> value;
	if (key.empty() || value.empty())
		throw std::runtime_error("Cannot parse input line: " + line);
	apply_key_value(config, current_test, key, value);
}

void parse_test_line(InputConfig& config, const std::string& line) {
	std::istringstream in(line);
	std::string marker;
	std::string name;
	in >> marker >> name;
	if (name.empty())
		throw std::runtime_error("Expected test name after '" + marker + "'.");

	TestRunConfig& run = append_test(config, name);
	std::string token;
	while (in >> token) {
		const size_t eq = token.find('=');
		if (eq == std::string::npos)
			throw std::runtime_error("Expected key=value in test line: " + token);
		apply_override(run, token.substr(0, eq), token.substr(eq + 1));
	}
}

InputConfig read_input(const std::filesystem::path& file_name) {
	InputConfig config;
	std::ifstream in(file_name);
	if (!in.is_open()) {
		std::cout << "input.txt was not found; all tests will run with default parameters.\n";
		add_test_name(config, "all");
		return config;
	}

	TestRunConfig* current_test = nullptr;
	std::string line;
	size_t line_number = 0;
	while (std::getline(in, line)) {
		++line_number;
		const size_t comment = line.find('#');
		if (comment != std::string::npos)
			line = line.substr(0, comment);

		line = trim(line);
		if (line.empty())
			continue;

		try {
			if (line.front() == '[' && line.back() == ']') {
				const std::string name = trim(line.substr(1, line.size() - 2));
				current_test = &append_test(config, name);
			}
			else if (starts_with(normalize(line), "test ") || starts_with(normalize(line), "тест ")) {
				parse_test_line(config, line);
				current_test = nullptr;
			}
			else {
				parse_assignment_line(config, current_test, line);
			}
		}
		catch (const std::exception& e) {
			throw std::runtime_error("input.txt:" + std::to_string(line_number) + ": " + e.what());
		}
	}

	if (config.tests.empty())
		add_test_name(config, "all");

	return config;
}

MeshSize apply_overrides(const TestCase& test, const TestRunConfig& run) {
	MeshSize mesh = test.mesh;
	if (run.bottom_x) mesh.bottom_x = *run.bottom_x;
	if (run.bottom_y) mesh.bottom_y = *run.bottom_y;
	if (run.top_x) mesh.top_x = *run.top_x;
	if (run.top_y) mesh.top_y = *run.top_y;
	if (run.lambda_nodes) mesh.lambda_nodes = *run.lambda_nodes;
	return mesh;
}

void validate_mesh(const MeshSize& mesh, ContactMethod contact_method) {
	if (mesh.bottom_x < 2 || mesh.bottom_y < 2 || mesh.top_x < 2 || mesh.top_y < 2)
		throw std::runtime_error("Every body mesh must have at least two nodes in each direction.");
	if (contact_method != ContactMethod::SlaveNodes && mesh.lambda_nodes < 2)
		throw std::runtime_error("lambda_nodes must be at least 2 for uniform contact methods.");
}

void validate_geometry(CoordinateSystem coordinate_system, const TestCase& test) {
	if (coordinate_system == CoordinateSystem::Axisymmetric &&
		(test.bottom_a.x <= 0.0 || test.top_a.x <= 0.0)) {
		throw std::runtime_error("Axisymmetric tests require positive radial coordinates.");
	}
}

std::string run_label(
	CoordinateSystem coordinate_system,
	ContactMethod contact_method,
	ContactSlaveBody slave_body,
	const MeshSize& mesh) {

	std::ostringstream label;
	label << "coord_" << to_string(coordinate_system)
		<< "__contact_" << to_string(contact_method);
	if (contact_method == ContactMethod::SlaveNodes)
		label << "__passive_" << to_string(slave_body);
	label << "__bottom_" << mesh.bottom_x << "x" << mesh.bottom_y
		<< "__top_" << mesh.top_x << "x" << mesh.top_y;
	if (contact_method != ContactMethod::SlaveNodes)
		label << "__lambda_" << mesh.lambda_nodes;
	return label.str();
}

double point_component(const Point& p, char component) {
	return (component == 'x' || component == 'r') ? p.x : p.y;
}

double derivative_x(
	const FEM& mesh,
	const std::vector<Point>& field,
	size_t row,
	size_t col,
	char component) {

	const size_t mx = mesh.xsize();
	auto value = [&](size_t c) {
		return point_component(field[row * mx + c], component);
	};

	if (mx == 2) {
		const double dx = mesh[row * mx + 1].x - mesh[row * mx].x;
		return (value(1) - value(0)) / dx;
	}

	if (col == 0) {
		const double dx = mesh[row * mx + 1].x - mesh[row * mx].x;
		return (-3.0 * value(0) + 4.0 * value(1) - value(2)) / (2.0 * dx);
	}

	if (col + 1 == mx) {
		const double dx = mesh[row * mx + col].x - mesh[row * mx + col - 1].x;
		return (3.0 * value(col) - 4.0 * value(col - 1) + value(col - 2)) / (2.0 * dx);
	}

	const double dx = mesh[row * mx + col + 1].x - mesh[row * mx + col - 1].x;
	return (value(col + 1) - value(col - 1)) / dx;
}

double derivative_y(
	const FEM& mesh,
	const std::vector<Point>& field,
	size_t row,
	size_t col,
	char component) {

	const size_t mx = mesh.xsize();
	const size_t ny = mesh.ysize();
	auto value = [&](size_t r) {
		return point_component(field[r * mx + col], component);
	};

	if (ny == 2) {
		const double dy = mesh[mx + col].y - mesh[col].y;
		return (value(1) - value(0)) / dy;
	}

	if (row == 0) {
		const double dy = mesh[mx + col].y - mesh[col].y;
		return (-3.0 * value(0) + 4.0 * value(1) - value(2)) / (2.0 * dy);
	}

	if (row + 1 == ny) {
		const double dy = mesh[row * mx + col].y - mesh[(row - 1) * mx + col].y;
		return (3.0 * value(row) - 4.0 * value(row - 1) + value(row - 2)) / (2.0 * dy);
	}

	const double dy = mesh[(row + 1) * mx + col].y - mesh[(row - 1) * mx + col].y;
	return (value(row + 1) - value(row - 1)) / dy;
}

double axisymmetric_hoop_strain(
	const FEM& mesh,
	const std::vector<Point>& field,
	size_t row,
	size_t col) {

	const size_t node_id = row * mesh.xsize() + col;
	const double r = mesh[node_id].x;
	if (std::fabs(r) < kTraceEps)
		return derivative_x(mesh, field, row, col, 'r');
	return field[node_id].x / r;
}

std::vector<double> side_axis_coordinates(const FSEM& body, char side, bool fem_nodes) {
	const auto nodes = fem_nodes ? body.get_side_fem_nodes(side) : body.get_side_nodes(side);
	std::vector<double> x(nodes.size());
	for (size_t i = 0; i < nodes.size(); ++i)
		x[i] = fem_nodes ? body.fem[nodes[i]].x : body[nodes[i]].x;
	return x;
}

size_t find_trace_segment(const std::vector<double>& x_nodes, double x) {
	if (x_nodes.size() < 2)
		throw std::runtime_error("At least two trace nodes are required.");

	if (x <= x_nodes.front() + kTraceEps)
		return 0;
	if (x >= x_nodes.back() - kTraceEps)
		return x_nodes.size() - 2;

	for (size_t i = 0; i + 1 < x_nodes.size(); ++i)
		if (x >= x_nodes[i] - kTraceEps && x <= x_nodes[i + 1] + kTraceEps)
			return i;

	throw std::runtime_error("Trace interpolation point is outside the contact interval.");
}

double interpolate_trace_value(
	const std::vector<double>& x_nodes,
	const std::vector<double>& values,
	double x) {

	const size_t segment = find_trace_segment(x_nodes, x);
	const double x_left = x_nodes[segment];
	const double x_right = x_nodes[segment + 1];
	const double value_left = values[segment];
	const double value_right = values[segment + 1];

	if (std::fabs(x_left - x_right) < kTraceEps)
		return value_left;

	const double t = (x - x_left) / (x_right - x_left);
	return (1.0 - t) * value_left + t * value_right;
}

std::vector<double> contact_output_grid(
	const std::vector<double>& bottom_x,
	const std::vector<double>& top_x) {

	const double contact_left = std::max(bottom_x.front(), top_x.front());
	const double contact_right = std::min(bottom_x.back(), top_x.back());

	std::vector<double> grid;
	grid.reserve(bottom_x.size() + top_x.size() + 2);
	grid.push_back(contact_left);
	grid.push_back(contact_right);

	auto append = [&](const std::vector<double>& nodes) {
		for (double x : nodes)
			if (x >= contact_left - kTraceEps && x <= contact_right + kTraceEps)
				grid.push_back(x);
	};

	append(bottom_x);
	append(top_x);

	std::sort(grid.begin(), grid.end());
	grid.erase(std::unique(grid.begin(), grid.end(),
		[](double lhs, double rhs) { return std::fabs(lhs - rhs) < kTraceEps; }),
		grid.end());
	return grid;
}

std::vector<double> recover_side_normal_stress(
	CoordinateSystem coordinate_system,
	const FSEM& body,
	const std::vector<Point>& field,
	char side,
	double E,
	double nu) {

	const auto side_nodes = body.get_side_fem_nodes(side);
	const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
	const double mu = E / (2.0 * (1.0 + nu));
	const size_t mx = body.fem.xsize();

	std::vector<double> sigma(side_nodes.size(), 0.0);
	for (size_t i = 0; i < side_nodes.size(); ++i) {
		const size_t node_id = side_nodes[i];
		const size_t row = node_id / mx;
		const size_t col = node_id % mx;

		if (coordinate_system == CoordinateSystem::Axisymmetric) {
			const double dur_dr = derivative_x(body.fem, field, row, col, 'r');
			const double duz_dz = derivative_y(body.fem, field, row, col, 'z');
			const double hoop_strain = axisymmetric_hoop_strain(body.fem, field, row, col);
			sigma[i] = lambda * (dur_dr + hoop_strain) + (lambda + 2.0 * mu) * duz_dz;
		}
		else {
			const double dux_dx = derivative_x(body.fem, field, row, col, 'x');
			const double duy_dy = derivative_y(body.fem, field, row, col, 'y');
			sigma[i] = lambda * dux_dx + (lambda + 2.0 * mu) * duy_dy;
		}
	}

	return sigma;
}

std::vector<double> recover_side_normal_displacement(
	const FSEM& body,
	const std::vector<Point>& field,
	char side) {

	const auto side_nodes = body.get_side_fem_nodes(side);
	std::vector<double> values(side_nodes.size(), 0.0);
	for (size_t i = 0; i < side_nodes.size(); ++i)
		values[i] = field[side_nodes[i]].y;
	return values;
}

void save_displacement_component(
	const std::filesystem::path& file_name,
	const FEM& mesh,
	const std::vector<Point>& field,
	char component) {

	std::ofstream out(file_name);
	out << std::setprecision(16);
	for (size_t i = 0; i < std::min(mesh.psize(), field.size()); ++i)
		out << mesh[i].x << " " << mesh[i].y << " "
			<< point_component(field[i], component) << "\n";
}

void save_displacement_vector(
	const std::filesystem::path& file_name,
	const FEM& mesh,
	const std::vector<Point>& field) {

	std::ofstream out(file_name);
	out << std::setprecision(16);
	for (size_t i = 0; i < std::min(mesh.psize(), field.size()); ++i)
		out << mesh[i].x << " " << mesh[i].y << " " << field[i].x << " " << field[i].y << "\n";
}

void save_contact_normal_stress(
	const std::filesystem::path& file_name,
	CoordinateSystem coordinate_system,
	const FSEM& bottom,
	const std::vector<Point>& bottom_field,
	const FSEM& top,
	const std::vector<Point>& top_field,
	double E,
	double nu) {

	const std::vector<double> bottom_x = side_axis_coordinates(bottom, 'N', true);
	const std::vector<double> top_x = side_axis_coordinates(top, 'S', true);
	const std::vector<double> grid = contact_output_grid(bottom_x, top_x);
	const std::vector<double> bottom_sigma =
		recover_side_normal_stress(coordinate_system, bottom, bottom_field, 'N', E, nu);
	const std::vector<double> top_sigma =
		recover_side_normal_stress(coordinate_system, top, top_field, 'S', E, nu);

	std::ofstream out(file_name);
	out << std::setprecision(16);
	for (double x : grid)
		out << x << " "
			<< interpolate_trace_value(bottom_x, bottom_sigma, x) << " "
			<< interpolate_trace_value(top_x, top_sigma, x) << "\n";
}

void save_contact_normal_displacement(
	const std::filesystem::path& file_name,
	const FSEM& bottom,
	const std::vector<Point>& bottom_field,
	const FSEM& top,
	const std::vector<Point>& top_field) {

	const std::vector<double> bottom_x = side_axis_coordinates(bottom, 'N', true);
	const std::vector<double> top_x = side_axis_coordinates(top, 'S', true);
	const std::vector<double> grid = contact_output_grid(bottom_x, top_x);
	const std::vector<double> bottom_u = recover_side_normal_displacement(bottom, bottom_field, 'N');
	const std::vector<double> top_u = recover_side_normal_displacement(top, top_field, 'S');

	std::ofstream out(file_name);
	out << std::setprecision(16);
	for (double x : grid) {
		const double ub = interpolate_trace_value(bottom_x, bottom_u, x);
		const double ut = interpolate_trace_value(top_x, top_u, x);
		out << x << " " << ub << " " << ut << " " << (ub - ut) << "\n";
	}
}

std::vector<double> lambda_nodes_for_output(
	const FSEM& bottom,
	const FSEM& top,
	const ContactOptions& options) {

	const std::vector<double> bottom_x = side_axis_coordinates(bottom, 'N', false);
	const std::vector<double> top_x = side_axis_coordinates(top, 'S', false);
	const double contact_left = std::max(bottom_x.front(), top_x.front());
	const double contact_right = std::min(bottom_x.back(), top_x.back());

	ContactDiscretization discretization;
	if (options.method == ContactMethod::SlaveNodes) {
		discretization = contact_slave_nodes::build(
			bottom_x, top_x, contact_left, contact_right, options.slave_body);
	}
	else if (options.method == ContactMethod::UniformLambdaPartition) {
		discretization = contact_uniform_lambda_partition::build(
			bottom_x, top_x, contact_left, contact_right, options.lambda_node_count);
	}
	else {
		discretization = contact_uniform_union_partition::build(
			bottom_x, top_x, contact_left, contact_right, options.lambda_node_count);
	}
	return discretization.lambda_nodes;
}

void save_lagrange_multipliers(
	const std::filesystem::path& file_name,
	const std::vector<double>& solution,
	size_t start,
	const std::vector<double>& lambda_nodes) {

	std::ofstream out(file_name);
	out << std::setprecision(16);
	const size_t count = std::min(lambda_nodes.size(), solution.size() - start);
	for (size_t i = 0; i < count; ++i)
		out << lambda_nodes[i] << " " << solution[start + i] << "\n";
}

double relative_error(
	const FEM& mesh,
	const std::vector<Point>& field,
	const vec_function& exact) {

	double numerator = 0.0;
	double denominator = 0.0;
	for (size_t i = 0; i < std::min(mesh.psize(), field.size()); ++i) {
		const Point u = exact(mesh[i]);
		const Point diff = field[i] - u;
		numerator += diff.x * diff.x + diff.y * diff.y;
		denominator += u.x * u.x + u.y * u.y;
	}

	return denominator > 1e-30 ? std::sqrt(numerator / denominator) : std::sqrt(numerator);
}

void save_parameters(
	const std::filesystem::path& file_name,
	const std::string& test_name,
	const MeshSize& mesh,
	const ContactOptions& options,
	double E,
	double nu) {

	std::ofstream out(file_name);
	out << "test=" << test_name << "\n";
	out << "coordinate=" << to_string(options.coordinate_system) << "\n";
	out << "contact=" << to_string(options.method) << "\n";
	out << "passive_body=" << to_string(options.slave_body) << "\n";
	out << "bottom_mesh=" << mesh.bottom_x << "x" << mesh.bottom_y << "\n";
	out << "top_mesh=" << mesh.top_x << "x" << mesh.top_y << "\n";
	out << "lambda_nodes=" << mesh.lambda_nodes << "\n";
	out << "E=" << std::setprecision(16) << E << "\n";
	out << "nu=" << nu << "\n";
}

void apply_dirichlet(FSEM& body, const std::array<char, 3>& sides, const vec_function& exact) {
	for (char side : sides)
		body.set_bc1(side, exact);
}

void run_test(const TestRunConfig& run, const InputConfig& config) {
	const auto it = TESTS.find(run.name);
	if (it == TESTS.end()) {
		std::cerr << "Unknown test '" << run.name << "'. It is skipped.\n";
		return;
	}

	const TestCase& test = it->second;
	const MeshSize mesh = apply_overrides(test, run);
	validate_mesh(mesh, config.contact_method);
	validate_geometry(config.coordinate_system, test);

	const vec_function exact = test.exact_solution(test.nu, test.E);
	FSEM bottom(test.E, test.nu, test.bottom_a, test.bottom_b,
		mesh.bottom_x, mesh.bottom_y, 1, 1, config.coordinate_system);
	FSEM top(test.E, test.nu, test.top_a, test.top_b,
		mesh.top_x, mesh.top_y, 1, 1, config.coordinate_system);

	bottom.construct_basis();
	top.construct_basis();
	apply_dirichlet(bottom, test.bottom_dirichlet_sides, exact);
	apply_dirichlet(top, test.top_dirichlet_sides, exact);

	ContactOptions contact_options;
	contact_options.coordinate_system = config.coordinate_system;
	contact_options.method = config.contact_method;
	contact_options.slave_body = config.slave_body;
	contact_options.lambda_node_count = mesh.lambda_nodes;

	const std::vector<double> solution = solve_mortar_contact(
		bottom, top, bottom.get_f(), top.get_f(), contact_options);

	const size_t n_bottom = bottom.get_K().size();
	const size_t n_top = top.get_K().size();
	const std::vector<Point> bottom_field = bottom.find_answer(solution);
	const std::vector<Point> top_field = top.find_answer(solution, n_bottom);

	const std::filesystem::path output_dir =
		std::filesystem::path("res") / run.name /
		run_label(config.coordinate_system, config.contact_method, config.slave_body, mesh);
	std::filesystem::create_directories(output_dir);

	const char first_component = config.coordinate_system == CoordinateSystem::Axisymmetric ? 'r' : 'x';
	const char second_component = config.coordinate_system == CoordinateSystem::Axisymmetric ? 'z' : 'y';

	save_displacement_vector(output_dir / "bottom_displacement.txt", bottom.fem, bottom_field);
	save_displacement_vector(output_dir / "top_displacement.txt", top.fem, top_field);
	save_displacement_component(output_dir / ("bottom_displacement_" + std::string(1, first_component) + ".txt"),
		bottom.fem, bottom_field, first_component);
	save_displacement_component(output_dir / ("bottom_displacement_" + std::string(1, second_component) + ".txt"),
		bottom.fem, bottom_field, second_component);
	save_displacement_component(output_dir / ("top_displacement_" + std::string(1, first_component) + ".txt"),
		top.fem, top_field, first_component);
	save_displacement_component(output_dir / ("top_displacement_" + std::string(1, second_component) + ".txt"),
		top.fem, top_field, second_component);

	save_contact_normal_stress(output_dir / "contact_normal_stress.txt",
		config.coordinate_system, bottom, bottom_field, top, top_field, test.E, test.nu);
	save_contact_normal_displacement(output_dir / "contact_normal_displacement.txt",
		bottom, bottom_field, top, top_field);
	save_lagrange_multipliers(output_dir / "lagrange_multipliers.txt",
		solution, n_bottom + n_top, lambda_nodes_for_output(bottom, top, contact_options));

	const double bottom_error = relative_error(bottom.fem, bottom_field, exact);
	const double top_error = relative_error(top.fem, top_field, exact);
	{
		std::ofstream out(output_dir / "error.txt");
		out << std::setprecision(16)
			<< "bottom_relative_error=" << bottom_error << "\n"
			<< "top_relative_error=" << top_error << "\n";
	}
	save_parameters(output_dir / "parameters.txt", run.name, mesh, contact_options, test.E, test.nu);

	std::cout << "Saved " << run.name << " -> " << output_dir.string()
		<< " (lambda=" << (solution.size() - n_bottom - n_top)
		<< ", error bottom=" << bottom_error
		<< ", top=" << top_error << ")\n";
}

} // namespace

int main() {
	try {
		const InputConfig config = read_input("input.txt");
		std::cout << "Coordinate system: " << to_string(config.coordinate_system) << "\n";
		std::cout << "Contact method: " << to_string(config.contact_method) << "\n";
		if (config.contact_method == ContactMethod::SlaveNodes)
			std::cout << "Passive body: " << to_string(config.slave_body) << "\n";

		for (const TestRunConfig& run : config.tests) {
			try {
				run_test(run, config);
			}
			catch (const std::exception& e) {
				std::cerr << "Test '" << run.name << "' failed: " << e.what() << "\n";
			}
		}
	}
	catch (const std::exception& e) {
		std::cerr << "Fatal error: " << e.what() << "\n";
		return 1;
	}

	return 0;
}
