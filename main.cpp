#include <iostream>
#include <fstream>
#include <vector>

#include "matrix.h"
#include "linear_equation.h"
#include "tests.h"

int main()
{
    Square_Matrix<double> A(20);
    for(size_t i = 0; i < 20; ++i)
        for(size_t j = 0; j < 20; ++j) {
            if(i == j)
                A(i, j) = 10.0;
            else
                A(i, j) = 1.0 / (i + j + 2.0);
        }
    std::vector<double> f(20);
    for(size_t i = 0; i < 20; ++i) {
        f[i] = 1.0 / (i + 1.0);
    }
    std::vector<double> x0(20);
    for(size_t i = 0; i < 20; ++i) {
        x0[i] = 1.0;
    }

    std::ofstream fout;
    fout.open("1.txt");

    auto sol_gauss = solve_linear_equation_by_Gauss(A, f);
    double gauss_nev = N1(A * sol_gauss - f);
    fout << "Gauss: (";
    for(const auto& el : sol_gauss) {
        fout << el << " ";
    }
    fout << ")" << std::endl;
    fout << "Gauss_nev: " << gauss_nev << std::endl;

    auto sol_zeidel = solve_linear_equation_by_Zeidel(A, f, x0);
    double zeidel_nev = N1(A * sol_zeidel - f);
    fout << "Zeidel: (";
    for(const auto& el : sol_zeidel) {
        fout << el << " ";
    }
    fout << ")" << std::endl;
    fout << "Zeidel_nev: " << zeidel_nev << std::endl;
    return 0;
}
