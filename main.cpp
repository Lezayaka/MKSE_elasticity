#include <iostream>
#include <numbers>
#include <cmath>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdexcept>

#define M_PI 3.14159265358979323846

#include "Headers.h"

vec_function get_func(double nu, double E) {
    (void)nu;
    (void)E;
    return [nu, E](const Point& p) {
        (void)nu;
        (void)E;
        (void)p;
        return Point{ 0, 0 };
        };
}

void save_displacement_component(
    const std::string& file_name,
    const FEM& mesh,
    const std::vector<Point>& displacement_field,
    char component) {

    std::ofstream out(file_name);
    if (!out.is_open()) {
        std::cerr << "Failed to open output file: " << file_name << '\n';
        return;
    }

    out << std::setprecision(16);
    const size_t n = std::min(mesh.psize(), displacement_field.size());
    for (size_t i = 0; i < n; ++i) {
        const Point& p = mesh[i];
        const Point& u = displacement_field[i];
        const bool radial_component = (component == 'x' || component == 'r');
        const double value = radial_component ? u.x : u.y;
        out << p.x << " " << p.y << " " << value << "\n";
    }
}

namespace {

constexpr double kStressTraceEps = 1e-12;
constexpr double kAxisEps = 1e-12;

bool almost_equal(double lhs, double rhs, double eps = kStressTraceEps) {
    return std::fabs(lhs - rhs) < eps;
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
    if (ny < 2)
        throw std::runtime_error("At least two grid nodes along y are required.");

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

    const size_t mx = mesh.xsize();
    const size_t node_id = row * mx + col;
    const double r = mesh[node_id].x;
    if (std::fabs(r) < kAxisEps)
        return derivative_x(mesh, field, row, col, 'r');

    return field[node_id].x / r;
}

std::vector<double> get_side_r_coordinates(const FSEM& body, char side) {
    const auto side_nodes = body.get_side_fem_nodes(side);
    std::vector<double> side_r(side_nodes.size());

    for (size_t i = 0; i < side_nodes.size(); ++i)
        side_r[i] = body.fem[side_nodes[i]].x;

    return side_r;
}

std::vector<double> recover_side_sigma_zz(
    const FSEM& body,
    const std::vector<Point>& displacement_field,
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

        const double dur_dr = derivative_x(body.fem, displacement_field, row, col, 'r');
        const double duz_dz = derivative_y(body.fem, displacement_field, row, col, 'z');
        const double hoop_strain = axisymmetric_hoop_strain(body.fem, displacement_field, row, col);
        sigma[i] = lambda * (dur_dr + hoop_strain) + (lambda + 2.0 * mu) * duz_dz;
    }

    return sigma;
}

std::vector<double> recover_side_normal_displacement(
    const FSEM& body,
    const std::vector<Point>& displacement_field,
    char side) {

    const auto side_nodes = body.get_side_fem_nodes(side);
    const char component =
        (side == 'N' || side == 'S') ? 'z' :
        (side == 'W' || side == 'E') ? 'r' : '\0';

    std::vector<double> values(side_nodes.size(), 0.0);
    for (size_t i = 0; i < side_nodes.size(); ++i)
        values[i] = point_component(displacement_field[side_nodes[i]], component);

    return values;
}

size_t find_trace_segment(const std::vector<double>& x_nodes, double x) {
    if (x_nodes.size() < 2)
        throw std::runtime_error("At least two trace nodes are required for interpolation.");

    if (x <= x_nodes.front() + kStressTraceEps)
        return 0;

    if (x >= x_nodes.back() - kStressTraceEps)
        return x_nodes.size() - 2;

    for (size_t i = 0; i + 1 < x_nodes.size(); ++i)
        if (x >= x_nodes[i] - kStressTraceEps && x <= x_nodes[i + 1] + kStressTraceEps)
            return i;

    throw std::runtime_error("Trace point is outside the interpolation range.");
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

    if (almost_equal(x_left, x_right))
        return value_left;

    const double t = (x - x_left) / (x_right - x_left);
    return (1.0 - t) * value_left + t * value_right;
}

std::vector<double> build_contact_x_grid(
    const std::vector<double>& bottom_x,
    const std::vector<double>& top_x) {
    if (bottom_x.empty() || top_x.empty())
        throw std::runtime_error("Empty contact boundary while preparing stress output.");

    const double contact_left = std::max(bottom_x.front(), top_x.front());
    const double contact_right = std::min(bottom_x.back(), top_x.back());
    if (contact_left > contact_right + kStressTraceEps)
        throw std::runtime_error("Bodies do not overlap along the contact trace.");

    std::vector<double> x_grid;
    x_grid.reserve(bottom_x.size() + top_x.size() + 2);
    x_grid.push_back(contact_left);
    x_grid.push_back(contact_right);

    auto append_contact_nodes = [&](const std::vector<double>& side_x) {
        for (double x : side_x)
            if (x >= contact_left - kStressTraceEps && x <= contact_right + kStressTraceEps)
                x_grid.push_back(x);
    };

    append_contact_nodes(bottom_x);
    append_contact_nodes(top_x);

    std::sort(x_grid.begin(), x_grid.end());
    x_grid.erase(
        std::unique(x_grid.begin(), x_grid.end(),
            [](double lhs, double rhs) { return almost_equal(lhs, rhs); }),
        x_grid.end());

    return x_grid;
}

void save_contact_normal_stress(
    const std::string& file_name,
    const FSEM& bottom_body,
    const std::vector<Point>& bottom_field,
    const FSEM& top_body,
    const std::vector<Point>& top_field,
    double E,
    double nu) {

    const std::vector<double> bottom_x = get_side_r_coordinates(bottom_body, 'N');
    const std::vector<double> top_x = get_side_r_coordinates(top_body, 'S');
    const std::vector<double> x_grid = build_contact_x_grid(bottom_x, top_x);

    const std::vector<double> sigma_bottom_nodes =
        recover_side_sigma_zz(bottom_body, bottom_field, 'N', E, nu);
    const std::vector<double> sigma_top_nodes =
        recover_side_sigma_zz(top_body, top_field, 'S', E, nu);

    std::ofstream out(file_name);
    if (!out.is_open()) {
        std::cerr << "Failed to open output file: " << file_name << '\n';
        return;
    }

    out << std::setprecision(16);
    for (double x : x_grid) {
        const double sigma_bottom = interpolate_trace_value(bottom_x, sigma_bottom_nodes, x);
        const double sigma_top = interpolate_trace_value(top_x, sigma_top_nodes, x);
        out << x << " " << sigma_bottom << " " << sigma_top << "\n";
    }
}

void save_contact_normal_displacement(
    const std::string& file_name,
    const FSEM& bottom_body,
    const std::vector<Point>& bottom_field,
    const FSEM& top_body,
    const std::vector<Point>& top_field) {

    const std::vector<double> bottom_x = get_side_r_coordinates(bottom_body, 'N');
    const std::vector<double> top_x = get_side_r_coordinates(top_body, 'S');
    const std::vector<double> x_grid = build_contact_x_grid(bottom_x, top_x);

    const std::vector<double> u_bottom_nodes =
        recover_side_normal_displacement(bottom_body, bottom_field, 'N');
    const std::vector<double> u_top_nodes =
        recover_side_normal_displacement(top_body, top_field, 'S');

    std::ofstream out(file_name);
    if (!out.is_open()) {
        std::cerr << "Failed to open output file: " << file_name << '\n';
        return;
    }

    out << std::setprecision(16);
    for (double x : x_grid) {
        const double u_bottom = interpolate_trace_value(bottom_x, u_bottom_nodes, x);
        const double u_top = interpolate_trace_value(top_x, u_top_nodes, x);
        out << x << " " << u_bottom << " " << u_top << "\n";
    }
}

}

//vec_function ans(double nu, double E) {
//    return [nu, E](const Point& p) {
//        double mu = E / (2 * (1 + nu));
//        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
//        const double& x = p.x, & y = p.y;
//
//        //     x2      y2      x        y        c
//        double a1 = 0, a3 = 0, a4 = 1, a5 = 0, a6 = 0;
//        double b1 = 0, b3 = 0, b4 = 0, b5 = 0, b6 = 0;
//
//        double a2 = -(2 * mu * b1 + (4 * mu + 2 * lambda) * b3) / (lambda + mu),
//            b2 = -(2 * mu * a3 + (4 * mu + 2 * lambda) * a1) / (lambda + mu);
//
//        return Point{
//            a1 * x * x + a2 * x * y + a3 * y * y + a4 * x + a5 * y + a6,
//            b1 * x * x + b2 * x * y + b3 * y * y + b4 * x + b5 * y + b6
//        };
//        };
//}

// u = {2r + 5 / r, 34z}
vec_function ans(double nu, double E) {
    return [nu, E](const Point& p) {
        (void)nu;
        (void)E;
        return Point{ 2 * p.x + 5 / p.x, 34 * p.y };
        };
}
/*vec_function ans(double nu, double E) {
    return [nu, E](const Point& p) {
        (void)nu;
        (void)E;
        return Point{ 2 * p.x + 5 / p.x, 34 * p.y};
        };
}*/

/*vec_function ans(double nu, double E) {
    return [nu, E](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        const double& x = p.x, & y = p.y;

        return Point{
            exp(x) * cos(y - 0.5),

            -exp(x)* sin(y - 0.5)
        };
        };
}*/

// пример для разных тел
/*vec_function ans(double nu, double E) {
    return [nu, E](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        const double& x = p.x, & y = p.y;

        return Point{
            -(lambda + 2 * mu) / (2 * lambda) * x * x -
            4 * (lambda + mu) / lambda * y +
            (3 * lambda + 4 * mu) / (2 * lambda) * y * y,

            x * y
        };
        };
}*/

// u = {3x^2 - 3 (y - 2)^2, -6x(y - 2)}
/*vec_function ans(double nu, double E) {
    return [nu, E](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        const double& x = p.x, & y = p.y;

        return Point{
            3 * x * x - 3 * (y - 2) * (y - 2),

            - 6 * x * (y - 2)
        };
        };
}*/

int main() {
    double E = 21e+10;
    double nu = 0.3;

    Point bottom_a = { 1, 0 }, bottom_b = { 3, 0.5 };
    Point top_a = { 1, 0.5 }, top_b = { 3, 3 };

    auto ANS = ans(nu, E);

    // несовпадающие сетки
    size_t col1 = 10, col2 = 3;
    size_t n_bottom_x = col1, n_bottom_y = col1, n_top_x = col1, n_top_y = col1;
    const size_t lambda_node_count = std::max(n_bottom_x, n_top_x);

    FSEM bottom(E, nu, bottom_a, bottom_b, n_bottom_x, n_bottom_y);
    FSEM top(E, nu, top_a, top_b, n_top_x, n_top_y);

    bottom.construct_basis();
    top.construct_basis();

    // Нижнее тело
    bottom.set_bc1('W', ANS); 

    bottom.set_bc1('E', ANS);
    bottom.set_bc1('S', ANS);

    // u = {x, 0}
    /*bottom.construct_f_bc2({ 1, 0, 0, 1 }, {
         [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{ - (lambda + 2 * mu), 0}; },

         zero,

        zero,

         [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{ 0, -lambda}; } // для (r, z) будет 2 * lambda
        });*/

    // u = {0, y}
    /*bottom.construct_f_bc2({ 0, 0, 1, 0 }, {
         zero,

         zero,

        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{ lambda, 0}; },

         zero
        });*/

    // u = {1 + 2x - 3y, 4 + 3x + 2y}
    /*bottom.construct_f_bc2({ 0, 0, 1, 0 }, {
        zero,
        zero,
        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{ 4 * (lambda + mu), 0}; },
        zero
        });*/

     // Верхнее тело
    top.set_bc1('W', ANS); 
    top.set_bc1('N', ANS);
    top.set_bc1('E', ANS);

    // u = {2r + 5 / r, 34z}
    /*top.construct_f_bc2({ 0, 1, 0, 0 }, {
         zero,

          [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{ 0, 2 * lambda * 2 + (lambda + 2 * mu) * 34}; },

        zero,

        zero
        });*/

    // u = {x, 0}
    /*top.construct_f_bc2({ 1, 1, 0, 0 }, {
         [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{ -(lambda + 2 * mu), 0}; },

          [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{ 0, lambda}; },

        zero,

        zero
        });*/

    // u = {0, y}
    /*top.construct_f_bc2({ 0, 0, 1, 0 }, {
         zero,

         zero,

        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{ lambda, 0}; },

         zero
        });*/

    // u = {-21x, 13y}
    /*top.construct_f_bc2({ 0, 1, 0, 0 }, {
         zero,

         [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{ 0, - 8 * lambda + 26 * mu }; },

        zero,

         zero
        });*/

    // u = {1 + 2x - 3y, 4 + 3x + 2y}
    /*top.construct_f_bc2({ 0, 0, 1, 0 }, {
        zero,
        zero,
        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{ 4 * (lambda + mu), 0}; },
        zero
        });*/

    // нагрузка для примера с разными телами, точным решением
   /*top.construct_f_bc2({ 0, 1, 0, 0 }, {
        zero,

        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        double y_bound = 3;
        return Point{ (y_bound - 1) * 4.0 * mu * (lambda + mu) / lambda, 0}; },

        zero,

        zero
        });*/

    // нагрузка для тестового примера
    /*top.construct_f_bc2({ 0, 1, 0, 0 }, {
        zero,

        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        double y_bound = 3;
        return Point{ 0, -2e+10 }; },

        zero,

        zero
        });*/

    std::vector<double> rhs_bottom = bottom.get_f();
    std::vector<double> rhs_top = top.get_f();

    std::vector<double> solution = solve_mortar_contact(
        bottom, top, rhs_bottom, rhs_top, lambda_node_count);

    const size_t n1 = bottom.get_K().size();
    const size_t n2 = top.get_K().size();

    std::cout << "Unknowns: u1=" << n1 / 2
        << " nodes, u2=" << n2 / 2
        << " nodes, lambda=" << (solution.size() - n1 - n2)
        << " nodes\n";

    auto res_bottom = bottom.find_answer(solution);
    auto res_top = top.find_answer(solution, n1);

    const auto bottom_contact_nodes = bottom.get_side_fem_nodes('N');
    const auto top_contact_nodes = top.get_side_fem_nodes('S');
    const size_t n_contact = std::min(bottom_contact_nodes.size(), top_contact_nodes.size());

    //std::cout << "\nContact surface values:\n";
    std::cout << "id\t(r, z)\tu_bottom\tu_top\n";

    for (size_t i = 0; i < n_contact; ++i) {
        const size_t bottom_contact_id = bottom_contact_nodes[i];
        const size_t top_contact_id = top_contact_nodes[i];

        const Point bottom_contact_point = bottom.fem[bottom_contact_id];
        const Point top_contact_point = top.fem[top_contact_id];
        //const Point exact_contact_value = ANS(bottom_contact_point);

        std::cout << i << "\t" << bottom_contact_point
            << "\t" << res_bottom[bottom_contact_id]
            << "\t" << res_top[top_contact_id]
            << "\n";
    }

    //std::cout << "\nu_bottom solution:\n";
    //for (size_t i = 0; i != res_bottom.size(); ++i) {
    //    //if (std::fabs(res_bottom[i].x - ANS((bottom.fem)[i]).x) > 1e-7 ||
    //      //  std::fabs(res_bottom[i].y - ANS((bottom.fem)[i]).y) > 1e-7)
    //    std::cout << i << ": " << res_bottom[i] << '\t' << ANS((bottom.fem)[i]) << '\n';
    //}

    //std::cout << "\nu_top solution:\n";
    //for (size_t i = 0; i != res_top.size(); ++i)
    //    //if (std::fabs(res_top[i].x - ANS((top.fem)[i]).x) > 1e-7 ||
    //      //  std::fabs(res_top[i].y - ANS((top.fem)[i]).y) > 1e-7)
    //        std::cout << i << ": " << res_top[i] << '\t' << ANS((top.fem)[i]) << '\n';

   save_displacement_component("results/bottom_displacement_r.txt", bottom.fem, res_bottom, 'r');
    save_displacement_component("results/bottom_displacement_z.txt", bottom.fem, res_bottom, 'z');
    save_displacement_component("results/top_displacement_r.txt", top.fem, res_top, 'r');
    save_displacement_component("results/top_displacement_z.txt", top.fem, res_top, 'z');
    save_contact_normal_stress(
        "results/contact_normal_stress.txt",
        bottom,
        res_bottom,
        top,
        res_top,
        E,
        nu);
    save_contact_normal_displacement(
        "results/contact_normal_displacement.txt",
        bottom,
        res_bottom,
        top,
        res_top);

    // Расчет и вывод нормы ошибки
    double bottom_numerator = 0.0, bottom_denominator = 0.0;

    for (size_t i = 0; i != res_bottom.size(); ++i) {
        Point U = ANS((bottom.fem)[i]);
        Point diff = res_bottom[i] - U;
        bottom_numerator += diff.x * diff.x + diff.y * diff.y;
        bottom_denominator += U.x * U.x + U.y * U.y;
    }
    double rel_bottom = (bottom_denominator > 1e-30)
        ? sqrt(bottom_numerator / bottom_denominator)
        : 0.0;

    std::cout << "Reletive u_bottom: " << rel_bottom << "\n";

    double top_numerator = 0.0, top_denominator = 0.0;

    for (size_t i = 0; i != res_top.size(); ++i) {
        Point U = ANS((top.fem)[i]);
        Point diff = res_top[i] - U;
        top_numerator += diff.x * diff.x + diff.y * diff.y;
        top_denominator += U.x * U.x + U.y * U.y;
    }
    double rel_top = (top_denominator > 1e-30)
        ? sqrt(top_numerator / top_denominator)
        : 0.0;

    std::cout << "Reletive u_top: " << rel_top << "\n";

    return 0;
}
