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
    return [nu, E](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        const double& x = p.x, & y = p.y;

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
        const double value = (component == 'x') ? u.x : u.y;
        out << p.x << " " << p.y << " " << value << "\n";
    }
}

void save_displacement_contact() {
}

namespace {

constexpr double kStressTraceEps = 1e-12;

bool almost_equal(double lhs, double rhs, double eps = kStressTraceEps) {
    return std::fabs(lhs - rhs) < eps;
}

double signed_double_area(const Point& p1, const Point& p2, const Point& p3) {
    return (p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y);
}

double point_component(const Point& p, char component) {
    return (component == 'x') ? p.x : p.y;
}

double derivative_x(
    const FEM& mesh,
    const std::vector<Point>& field,
    size_t row,
    size_t col,
    char component) {

    const size_t mx = mesh.xsize();
    if (mx < 2)
        throw std::runtime_error("At least two grid nodes along x are required.");

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

std::vector<double> get_side_x_coordinates(const FSEM& body, char side) {
    const auto side_nodes = body.get_side_fem_nodes(side);
    std::vector<double> side_x(side_nodes.size());

    for (size_t i = 0; i < side_nodes.size(); ++i)
        side_x[i] = body.fem[side_nodes[i]].x;

    return side_x;
}

std::vector<double> recover_side_sigma_yy(
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

        const double dux_dx = derivative_x(body.fem, displacement_field, row, col, 'x');
        const double duy_dy = derivative_y(body.fem, displacement_field, row, col, 'y');
        sigma[i] = lambda * dux_dx + (lambda + 2.0 * mu) * duy_dy;
    }

    return sigma;
}

size_t find_trace_segment(const std::vector<double>& x_nodes, double x) {
    if (x_nodes.size() < 2)
        throw std::runtime_error("Contact trace has less than two nodes.");

    if (x <= x_nodes.front() + kStressTraceEps)
        return 0;

    if (x >= x_nodes.back() - kStressTraceEps)
        return x_nodes.size() - 2;

    for (size_t i = 0; i + 1 < x_nodes.size(); ++i)
        if (x >= x_nodes[i] - kStressTraceEps && x <= x_nodes[i + 1] + kStressTraceEps)
            return i;

    throw std::runtime_error("Requested x is outside the contact trace.");
}

double interpolate_trace_value(
    const std::vector<double>& x_nodes,
    const std::vector<double>& values,
    double x) {

    if (x_nodes.size() != values.size())
        throw std::runtime_error("Contact trace interpolation data is inconsistent.");

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
    if (contact_right - contact_left <= kStressTraceEps)
        throw std::runtime_error("Contact overlap is empty while preparing stress output.");

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

    const std::vector<double> bottom_x = get_side_x_coordinates(bottom_body, 'N');
    const std::vector<double> top_x = get_side_x_coordinates(top_body, 'S');
    const std::vector<double> x_grid = build_contact_x_grid(bottom_x, top_x);

    const std::vector<double> sigma_bottom_nodes =
        recover_side_sigma_yy(bottom_body, bottom_field, 'N', E, nu);
    const std::vector<double> sigma_top_nodes =
        recover_side_sigma_yy(top_body, top_field, 'S', E, nu);

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

}

// простые решения
/*vec_function ans(double nu, double E) {
    return [nu, E](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        const double& x = p.x, & y = p.y;

        //     x2      y2      x        y        c
        double a1 = 0, a3 = 0, a4 = 1, a5 = 0, a6 = 0;
        double b1 = 0, b3 = 0, b4 = 0, b5 = 0, b6 = 0;

        double a2 = -(2 * mu * b1 + (4 * mu + 2 * lambda) * b3) / (lambda + mu),
            b2 = -(2 * mu * a3 + (4 * mu + 2 * lambda) * a1) / (lambda + mu);

        return Point{
            a1 * x * x + a2 * x * y + a3 * y * y + a4 * x + a5 * y + a6,
            b1 * x * x + b2 * x * y + b3 * y * y + b4 * x + b5 * y + b6
        };
        };
}*/

// u = {e(x) cos(y - 0.5), - e(x) sin(y - 0.5)} 
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

    Point bottom_a = { 0, 0 }, bottom_b = { 3, 1 };
    Point top_a = { 0, 1 }, top_b = { 2, 4 };

    //auto ANS = ans(nu, E);

    // совпадающие сетки
    //size_t n_bottom_x = 3, n_bottom_y = 3, n_top_x = 3, n_top_y = 3;

    // несовпадающие сетки
    size_t n_bottom_x = 5, n_bottom_y = 5, n_top_x = 5, n_top_y = 5;
    //const bool bottom_is_master = false;
    const bool bottom_is_master = true;

    FSEM bottom(E, nu, bottom_a, bottom_b, n_bottom_x, n_bottom_y);
    FSEM top(E, nu, top_a, top_b, n_top_x, n_top_y);

    bottom.construct_basis();
    top.construct_basis();

    // Нижнее тело
    bottom.set_bc1('W', zero); 

    //bottom.set_bc1('E', zero);
    bottom.set_bc1('S', zero);

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
        return Point{ 0, -lambda}; }
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
    top.set_bc1('W', zero); 
    //top.set_bc1('N', zero);
    //top.set_bc1('E', zero);

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
    top.construct_f_bc2({ 0, 1, 0, 0 }, {
        zero,

        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        double y_bound = 3;
        return Point{ 0, -2e+10 }; },

        zero,

        zero
        });
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
        bottom, top, rhs_bottom, rhs_top, bottom_is_master);

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

    std::cout << "\nContact surface values:\n";
    std::cout << "id\tpoint\tu_bottom\tu_top\tu_exact(one-body)\n";

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

    save_displacement_component("results/bottom_displacement_x.txt", bottom.fem, res_bottom, 'x');
    save_displacement_component("results/bottom_displacement_y.txt", bottom.fem, res_bottom, 'y');
    save_displacement_component("results/top_displacement_x.txt", top.fem, res_top, 'x');
    save_displacement_component("results/top_displacement_y.txt", top.fem, res_top, 'y');
    save_contact_normal_stress(
        "results/contact_normal_stress.txt",
        bottom,
        res_bottom,
        top,
        res_top,
        E,
        nu);
    save_displacement_contact();

    // Расчет и вывод нормы ошибки
    /*double bottom_numerator = 0.0, bottom_denominator = 0.0;

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

    std::cout << "Reletive u_top: " << rel_top << "\n";*/

    return 0;
}
