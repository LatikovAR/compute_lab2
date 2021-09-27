#include "linear_equation.h"

std::vector<double> solve_linear_equation_by_Gauss(Matrix<double> A, std::vector<double> f)
{
    if(A.column_size() != f.size())
        throw std::invalid_argument("f non-suitable to A");
    if(A.row_size() > f.size())
        throw std::invalid_argument("probability of multiple solution of the equation");

    std::vector<size_t> vars_ordering(A.row_size());
    for(size_t i = 0; i < vars_ordering.size(); ++i)
        vars_ordering[i] = i;

    //--------------------------------straight pass----------------------------------------

    size_t column_size = A.column_size();
    for(size_t row_num = 0; row_num < column_size; ++row_num)
    {
        //choosing main element
        const auto& row = A.row(row_num);
        auto max_it = std::max_element(row.begin() + row_num, row.end());
        auto min_it = std::min_element(row.begin() + row_num, row.end());
        auto main_it = max_it;
        if(fabs(*max_it) < fabs(*min_it))
            main_it = min_it;
        if(main_it != row.begin() + row_num) {
            size_t max_pos = std::distance(row.begin(), max_it);
            A.swap_columns(row_num, max_pos);
            std::swap(vars_ordering[row_num], vars_ordering[max_pos]);
        }

        if(fabs(row[row_num]) < 1e-10) {
            if(fabs(f[row_num]) < 1e-10)
                throw std::invalid_argument("probability of multiple solution of the equation");
            else
                throw std::invalid_argument("equation doesn't have solution");
        }

        for(size_t i = row_num + 1; i != A.column_size(); ++i) {
            double d = A(i, row_num) / row[row_num];
            A(i, row_num) = 0.0;
            for(size_t j = row_num + 1; j < A.row_size(); ++j) {
                A(i, j) -= d * row[j];
            }
            f[i] -= d * f[row_num];
        }
    }

    //---------------------------------reverse pass----------------------------------------

    std::vector<double> solution(A.row_size());
    for(size_t i = A.row_size() - 1;; --i) {
        double& el = solution[i];
        el = f[i];
        for(size_t j = A.row_size() - 1; j > i; --j) {
            el -= solution[j] * A(i, j);
        }
        el /= A(i, i);

        if(i == 0) break;
    }

    //fix ordering
    std::vector<std::pair<size_t, double>> sort_solution(solution.size());
    for(size_t i = 0; i != solution.size(); ++i) {
        sort_solution[i].first = vars_ordering[i];
        sort_solution[i].second = solution[i];
    }
    std::sort(sort_solution.begin(), sort_solution.end(),
              [](const auto& lhs, const auto& rhs){
        return lhs.first < rhs.first;
    });

    std::vector<double> ordered_solution(sort_solution.size());
    for(size_t i = 0; i < sort_solution.size(); ++i)
        ordered_solution[i] = sort_solution[i].second;

    return ordered_solution;
}

std::vector<double> solve_linear_equation_by_Zeidel(const Square_Matrix<double>& A,
                                                    const std::vector<double>& f,
                                                    std::vector<double> x,
                                                    double accuracy)
{
    if(A.size() != f.size())
        throw std::invalid_argument("f non-suitable to A");
    if(f.size() != x.size())
        throw std::invalid_argument("x0 has wrong size");

    Lower_Triangle_Matrix<double> L(A.size());
    for(size_t i = 0; i < A.size(); ++i)
        for(size_t j = 0; j <= i; ++j)
            L(i, j) = A(i, j);
    Square_Matrix<double> U(A.size());
    for(size_t i = 0; i < A.size(); ++i)
        for(size_t j = i + 1; j < A.size(); ++j)
            U(i, j) = A(i, j);

    Square_Matrix<double> I = Inverse(L);
    Square_Matrix<double> B = -1.0 * I * U;
    std::vector<double> F = I * f;

    while(N1(A * x - f) > accuracy) {
        x = B * x + F;
    }

    return x;
}
