#include <iostream>
#include <numbers>
#include <cmath>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <string>

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

//vec_function ans(double nu, double E) {
//    return [nu, E](const Point& p) {
//        double mu = E / (2 * (1 + nu));
//        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
//        const double& x = p.x, & y = p.y;
//
//        return Point{
//            -((lambda + 2.0 * mu) / (2.0 * lambda)) * x * x
//           - (4.0 * (lambda + mu) / lambda) * y
//           + ((3.0 * lambda + 4.0 * mu) / (2.0 * lambda)) * y * y,
//
//            x * y
//        };
//        };
//}


int main() {
    double E = 21e+10;
    double nu = 0.3;

    Point bottom_a = { 0, 0 }, bottom_b = { 3, 1 };
    Point top_a = { 0, 1 }, top_b = { 2, 4 };

    //auto ANS = ans(nu, E);

    size_t n_x_b = 13, n_x_t = 9, n_y = 3;

    FSEM bottom(E, nu, bottom_a, bottom_b, n_x_b, 4);
    FSEM top(E, nu, top_a, top_b, n_x_t, 10);

    bottom.construct_basis();
    top.construct_basis();

    // Нижнее тело: фиксируем низ, остальные стороны свободны.
    bottom.set_bc1('W', [&](const Point& p) {
        return Point{ 0, NAN }; });

    //bottom.set_bc1('E', ans(nu, E));
    bottom.set_bc1('S', zero);

    /*bottom.construct_f_bc2({ 0, 0, 1, 0 }, {
        zero,

         zero,

         [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{ lambda * 3 + 2 * mu, 0}; },

         zero
     });*/

     // Верхнее тело: задаем внешнюю нагрузку сверху.
    top.set_bc1('W', [&](const Point& p) {
        return Point{ 0, NAN }; });
    // top.set_bc1('N', ans(nu, E));
    //top.set_bc1('E', ans(nu, E));


    top.construct_f_bc2({ 0, 1, 0, 0 }, {
        zero,

        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{ 0, -2e+10}; },

        zero,

        zero
        });

    std::vector<double> rhs_bottom = bottom.get_f();
    std::vector<double> rhs_top = top.get_f();

    std::vector<double> solution = solve_mortar_contact(bottom, top, rhs_bottom, rhs_top);

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

    /*std::cout << "\nContact surface values:\n";
    std::cout << "id\tpoint\tu_bottom\tu_top\tu_exact(one-body)\n";

    for (size_t i = 0; i < n_contact; ++i) {
        const size_t bottom_contact_id = bottom_contact_nodes[i];
        const size_t top_contact_id = top_contact_nodes[i];

        const Point bottom_contact_point = bottom.fem[bottom_contact_id];
        const Point top_contact_point = top.fem[top_contact_id];
        const Point exact_contact_value = ANS(bottom_contact_point);

        std::cout << i << "\t" << bottom_contact_point
            << "\t" << res_bottom[bottom_contact_id]
            << "\t" << res_top[top_contact_id]
            << "\t" << exact_contact_value << "\n";*/
    //}

    std::cout << "\nu_bottom solution:\n";
    for (size_t i = 0; i != res_bottom.size(); ++i) {
        //if (std::fabs(res_bottom[i].x - ANS((bottom.fem)[i]).x) > 1e-7 ||
          //  std::fabs(res_bottom[i].y - ANS((bottom.fem)[i]).y) > 1e-7)
        std::cout << i << ": " << res_bottom[i] << '\n';// << '\t' << ANS((bottom.fem)[i]) << '\n';
    }

    std::cout << "\nu_top solution:\n";
    for (size_t i = 0; i != res_top.size(); ++i)
        //if (std::fabs(res_top[i].x - ANS((top.fem)[i]).x) > 1e-7 ||
          //  std::fabs(res_top[i].y - ANS((top.fem)[i]).y) > 1e-7)
            std::cout << i << ": " << res_top[i] << '\n';// << '\t' << ANS((top.fem)[i]) << '\n';

    save_displacement_component("bottom_displacement_x.txt", bottom.fem, res_bottom, 'x');
    save_displacement_component("bottom_displacement_y.txt", bottom.fem, res_bottom, 'y');
    save_displacement_component("top_displacement_x.txt", top.fem, res_top, 'x');
    save_displacement_component("top_displacement_y.txt", top.fem, res_top, 'y');

    // Расчет и вывод нормы ошибки
   /* double bottom_numerator = 0.0, bottom_denominator = 0.0;

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