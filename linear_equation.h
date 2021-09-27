#pragma once

#include <stdexcept>
#include <algorithm>
#include <cmath>

#include "matrix.h"

std::vector<double> solve_linear_equation_by_Gauss(Matrix<double> A, std::vector<double> f);
std::vector<double> solve_linear_equation_by_Zeidel(const Square_Matrix<double> &A,
                                                    const std::vector<double> &f,
                                                    std::vector<double> x,
                                                    double accuracy = 1e-5);
