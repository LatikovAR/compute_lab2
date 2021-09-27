#include <iostream>

#include "tests.h"
#include "matrix.h"
#include "linear_equation.h"

void gauss_test() {
    std::cout << "Gauss test:\n";
    {
        Matrix A(2, 2);
        A(0, 0) = 1.0; A(0, 1) = 1.0;
        A(1, 0) = 1.0; A(1, 1) = -1.0;
        std::vector<double> f(2);
        f[0] = 1.0; f[1] = 0.0;
        auto sol = solve_linear_equation_by_Gauss(A, f);
        std::cout << "1: ";
        for(const auto& el : sol) {
            std::cout << el << " ";
        }
        std::cout << std::endl;
    }

    {
        Matrix A(4, 4);
        A(0, 0) = 3.0; A(0, 1) = 4.0; A(0, 2) = 2.0; A(0, 3) = 1.0;
        A(1, 0) = 1.0; A(1, 1) = 7.0; A(1, 2) = 1.0; A(1, 3) = 1.0;
        A(2, 0) = 2.0; A(2, 1) = 1.0; A(2, 2) = 3.0; A(2, 3) = 5.0;
        A(3, 0) = 4.0; A(3, 1) = -3.0; A(3, 2) = 4.0; A(3, 3) = 6.0;
        std::vector<double> f(4);
        f[0] = 16.0; f[1] = 23.0; f[2] = 10.0; f[3] = 1.0;
        auto sol = solve_linear_equation_by_Gauss(A, f);
        std::cout << "2: ";
        for(const auto& el : sol) {
            std::cout << el << " ";
        }
        std::cout << std::endl;
    }
}

void zeidel_test() {
    std::cout << "Zeidel test:\n";
    {
        Square_Matrix A(4);
        A(0, 0) = 3.0; A(0, 1) = 4.0; A(0, 2) = 2.0; A(0, 3) = 1.0;
        A(1, 0) = 1.0; A(1, 1) = 7.0; A(1, 2) = 1.0; A(1, 3) = 1.0;
        A(2, 0) = 2.0; A(2, 1) = 1.0; A(2, 2) = 3.0; A(2, 3) = 5.0;
        A(3, 0) = 4.0; A(3, 1) = -3.0; A(3, 2) = 4.0; A(3, 3) = 6.0;
        std::vector<double> f(4);
        f[0] = 16.0; f[1] = 23.0; f[2] = 10.0; f[3] = 1.0;
        std::vector<double> x0(4);
        x0[0] = 1.0; x0[1] = 1.0; x0[2] = 1.0; x0[3] = 1.0;
        auto sol = solve_linear_equation_by_Zeidel(A, f, x0);
        std::cout << "1: ";
        for(const auto& el : sol) {
            std::cout << el << " ";
        }
        std::cout << std::endl;
    }
}
