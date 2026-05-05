#include <cmath>
#include <iomanip>
#include <iostream>

#include "Headers.h"

vec_function ans(double nu, double E) {
    return [nu, E](const Point& p) {
        (void)nu;
        (void)E;
        return Point{ 2 * p.x + 5 / p.x, 34 * p.y };
    };
}

double relative_error(
    const FSEM& body,
    const std::vector<Point>& numerical,
    const vec_function& exact) {

    double numerator = 0.0;
    double denominator = 0.0;

    for (size_t i = 0; i < numerical.size(); ++i) {
        const Point exact_value = exact(body.fem[i]);
        const Point diff = numerical[i] - exact_value;
        numerator += diff.x * diff.x + diff.y * diff.y;
        denominator += exact_value.x * exact_value.x + exact_value.y * exact_value.y;
    }

    return denominator > 1e-30 ? std::sqrt(numerator / denominator) : 0.0;
}

int main() {
    const double E = 21e+10;
    const double nu = 0.3;

    const Point bottom_a = { 1, 0 };
    const Point bottom_b = { 3, 0.5 };
    const Point top_a = { 1, 0.5 };
    const Point top_b = { 3, 3 };

    const auto exact = ans(nu, E);

    std::cout << std::setprecision(16);
    std::cout << "col1 rel_bottom rel_top\n";

    for (size_t col1 = 3; col1 <= 20; ++col1) {
        const size_t n_bottom_x = col1;
        const size_t n_bottom_y = col1;
        const size_t n_top_x = col1;
        const size_t n_top_y = col1;
        const size_t lambda_node_count = col1;

        FSEM bottom(E, nu, bottom_a, bottom_b, n_bottom_x, n_bottom_y);
        FSEM top(E, nu, top_a, top_b, n_top_x, n_top_y);

        bottom.construct_basis();
        top.construct_basis();

        bottom.set_bc1('W', exact);
        bottom.set_bc1('E', exact);
        bottom.set_bc1('S', exact);

        top.set_bc1('W', exact);
        top.set_bc1('N', exact);
        top.set_bc1('E', exact);

        const std::vector<double> rhs_bottom = bottom.get_f();
        const std::vector<double> rhs_top = top.get_f();

        const std::vector<double> solution = solve_mortar_contact(
            bottom, top, rhs_bottom, rhs_top, lambda_node_count);

        const size_t n1 = bottom.get_K().size();
        const auto res_bottom = bottom.find_answer(solution);
        const auto res_top = top.find_answer(solution, n1);

        const double rel_bottom = relative_error(bottom, res_bottom, exact);
        const double rel_top = relative_error(top, res_top, exact);

        std::cout << col1 << " " << rel_bottom << " " << rel_top << "\n";
    }

    return 0;
}
