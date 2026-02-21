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
        double a1 = 0, a3 = 0, a4 = 10, a5 = 0,  a6 = 0;
        double b1 = 0, b3 = 0, b4 = 0,  b5 = -5, b6 = 0;

        double a2 = -(2 * mu * b1 + (4 * mu + 2 * lambda) * b3) / (lambda + mu),
            b2 = -(2 * mu * a3 + (4 * mu + 2 * lambda) * a1) / (lambda + mu);

        return Point{ 
            a1 * x * x + a2 * x * y + a3 * y * y + a4 * x + a5 * y + a6,
            b1* x* x + b2 * x * y + b3 * y * y + b4 * x + b5 * y + b6
        };
    };
}


int main()
{
    // Задаем условия задачи
    double E = 21e+10; // модуль Юнга
    double nu = 0.3; // коэффициент Пуассона
    Point a = { 0, 0 }, b = { 4, 2 }; // границы сетки
    size_t n_x = 4, n_y = 3; // количество узлов на границах области
    auto ANS = ans(nu, E);

    FSEM fsem(E, nu, a, b, n_x, n_y); // строим сетку на границе области
    
    //std::cout << "Nodes:\n";
    //fsem.print_nodes();
    fsem.construct_basis();
    
    fsem.set_bc1('W', ans(nu, E));
    //fsem.set_bc1('N', ans(nu, E));
    //fsem.set_bc1('E', ans(nu, E));
    fsem.set_bc1('S', ans(nu, E));

    fsem.set_bc2({ 0, 1, 1, 0 }, {
        zero, 

        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{0, lambda * 5 - 10 * mu}; },

        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{lambda * 5 + 2 * mu * 10, 0}; },

        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{-2 * 1.23 * mu, 0}; }
        });

    fsem.construct_f_bc2({ 0, 1, 1, 0 }, {
        zero,

        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{0, lambda * 5 - 10 * mu}; },

        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{lambda * 5 + 2 * mu * 10, 0}; },

        [&](const Point& p) {
        double mu = E / (2 * (1 + nu));
        double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        return Point{-2 * 1.23 * mu, 0}; }
        });

    auto res = fsem.find_answer();

   for (size_t i = 0; i != res.size(); ++i) {
        std::cout << res[i] << '\t' << ANS((fsem.fem)[i]) << '\n';
    }

    // Расчет и вывод нормы ошибки
    double max_x = 0, max_y = 0;
    
    for (size_t i = 0; i != res.size(); ++i) {
        Point U = ANS((fsem.fem)[i]);
        if (fabs(U.x) > 1e-15 and max_x < abs((U.x - res[i].x) / U.x))
            max_x = abs((U.x - res[i].x) / U.x);
        if (fabs(U.y) > 1e-15 and max_y < abs((U.y - res[i].y) / U.y))
            max_y = abs((U.y - res[i].y) / U.y);
    }

    std::cout << "Reletive: ";
    std::cout << std::max(max_x, max_y) << "\n";
}