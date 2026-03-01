#include <iostream>
#include <numbers>
#include <cmath>
#include <math.h>

# define M_PI 3.14159265358979323846

#include "Headers.h"

vec_function get_func(double nu, double E) {
    return [nu, E](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        const double& x = p.x, & y = p.y;

        return Point{0, 0};
    };
}

// Точное решение
//Point ans(const Point& p) {
//    return{ sin(M_PI*p.x), 0 };
//}

vec_function ans(double nu, double E) {
    return [nu, E](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        const double& x = p.x, & y = p.y;

        //     x2      y2      x        y        c
        double a1 = 0, a3 = 0, a4 = 0, a5 = 1,  a6 = 234;
        double b1 = 0, b3 = 0, b4 = 0,  b5 = 0, b6 = 5;

        double a2 = -(2 * mu * b1 + (4 * mu + 2 * lambda) * b3) / (lambda + mu),
            b2 = -(2 * mu * a3 + (4 * mu + 2 * lambda) * a1) / (lambda + mu);

        return Point{ 
            a1 * x * x + a2 * x * y + a3 * y * y + a4 * x + a5 * y + a6,
            b1* x* x + b2 * x * y + b3 * y * y + b4 * x + b5 * y + b6
        };
    };
}


int main() {
    double E = 21e+10;
    double nu = 0.3;

    Point bottom_a = { 0, 0 }, bottom_b = { 4, 2 };
    Point top_a = { 0, 2 }, top_b = { 4, 4 };

    size_t n_x = 5, n_y = 5;

    FSEM bottom(E, nu, bottom_a, bottom_b, n_x, n_y);
    FSEM top(E, nu, top_a, top_b, n_x, n_y);

    bottom.construct_basis();
    top.construct_basis();

    // Нижнее тело: фиксируем низ, остальные стороны свободны.
    bottom.set_bc1('S', ans(nu, E));
    bottom.set_bc1('W', ans(nu, E));
    bottom.set_bc1('E', ans(nu, E));

    // Верхнее тело: задаем внешнюю нагрузку сверху.
    top.set_bc1('W', ans(nu, E));
    top.set_bc1('N', ans(nu, E));
    top.set_bc1('E', ans(nu, E));

    std::vector<double> rhs_bottom = bottom.get_f();
    std::vector<double> rhs_top = top.get_f();
    
   std::vector<double> solution = solve_mortar_contact(bottom, top, rhs_bottom, rhs_top);
   
    const size_t n1 = bottom.get_K().size();
    const size_t n2 = top.get_K().size();
    
    std::cout << "Unknowns: u1=" << n1 / 2
        << " nodes, u2=" << n2 / 2
        << " nodes, lambda=" << (solution.size() - n1 - n2)
        << " nodes\n";

    auto ANS = ans(nu, E);

    auto res_bottom = bottom.find_answer(solution);
    auto res_top = top.find_answer(solution, n1);

    std::cout << "\nu_bottom solution:\n";
    for (size_t i = 0; i != res_bottom.size(); ++i)
        std::cout << i << ": " << res_bottom[i] << '\t' << ANS((bottom.fem)[i]) << '\n';

    std::cout << "\nu_top solution:\n";
    for (size_t i = 0; i != res_top.size(); ++i)
        std::cout << i << ": " << res_top[i] << '\t' << ANS((top.fem)[i]) << '\n';

    // Расчет и вывод нормы ошибки
    double max_x = 0, max_y = 0;

    for (size_t i = 0; i != res_bottom.size(); ++i) {
        Point U = ANS((bottom.fem)[i]);
        if (fabs(U.x) > 1e-15 and max_x < abs((U.x - res_bottom[i].x) / U.x))
            max_x = abs((U.x - res_bottom[i].x) / U.x);
        if (fabs(U.y) > 1e-15 and max_y < abs((U.y - res_bottom[i].y) / U.y))
            max_y = abs((U.y - res_bottom[i].y) / U.y);
    }

    std::cout << "Reletive u_bottom: ";
    std::cout << std::max(max_x, max_y) << "\n";

    max_x = 0, max_y = 0;

    for (size_t i = 0; i != res_top.size(); ++i) {
        Point U = ANS((top.fem)[i]);
        if (fabs(U.x) > 1e-15 and max_x < abs((U.x - res_top[i].x) / U.x))
            max_x = abs((U.x - res_top[i].x) / U.x);
        if (fabs(U.y) > 1e-15 and max_y < abs((U.y - res_top[i].y) / U.y))
            max_y = abs((U.y - res_top[i].y) / U.y);
    }

    std::cout << "Reletive u_top: ";
    std::cout << std::max(max_x, max_y) << "\n";

    return 0;
}