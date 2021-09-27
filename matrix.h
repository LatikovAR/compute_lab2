#pragma once

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>

template<typename T>
class Lower_Triangle_Matrix;

template<typename T = double>
class Matrix {
protected:
    std::vector<std::vector<T>> data_;

public:
    Matrix(size_t column_size, size_t row_size) :
        data_(column_size, std::vector<T>(row_size, static_cast<T>(0.0))) {}
    virtual ~Matrix() = default;

    Matrix(const Lower_Triangle_Matrix<T>& L) :
        data_(L.size(), std::vector<T>(L.size(), static_cast<T>(0.0))) {
        for(size_t i = 0; i < L.size(); ++i)
            for(size_t j = 0; j <= i; ++j)
                data_[i][j] = L(i, j);
    }

    T& operator()(size_t i, size_t j) {
        return data_.at(i).at(j);
    }
    const T& operator()(size_t i, size_t j) const {
        return data_.at(i).at(j);
    }

    size_t row_size() const {
        if(data_.size() == 0) return 0;
        return data_[0].size();
    }
    size_t column_size() const {
        return data_.size();
    }

    const std::vector<T>& row(size_t num) const {
        return data_.at(num);
    }

    void swap_columns(size_t col1_num, size_t col2_num) {
        std::for_each(data_.begin(), data_.end(), [=](auto& row){
            std::swap(row.at(col1_num), row.at(col2_num));
        });
    }

    void add_row_with_k(size_t src_num, size_t dst_num, T k = static_cast<T>(1.0)) {
        if(src_num >= column_size() || dst_num >= column_size())
            throw std::invalid_argument("Matrix: no that row");

        for(size_t i = 0; i < row_size(); ++i)
            data_[dst_num][i] += (k * data_[src_num][i]);
    }

    void mult_row_to_k(size_t row_num, T k) {
        if(row_num >= column_size())
            throw std::invalid_argument("Matrix: no that row");

        for(auto& el : data_[row_num]) {
            el *= k;
        }
    }

    Matrix<T>& operator*=(const Matrix<T>& rhs) {
        if(row_size() != rhs.column_size())
            throw std::invalid_argument("Matrix*: wrong sizes");
        Matrix<T> A_new(column_size(), rhs.row_size());

        for(size_t i = 0; i < A_new.column_size(); ++i)
            for(size_t j = 0; j < A_new.row_size(); ++j)
                for(size_t k = 0; k < row_size(); ++k)
                    A_new(i, j) += (data_[i][k] * rhs(k, j));

        data_ = std::move(A_new.data_);
        return *this;
    }
};

template<typename T = double>
Matrix<T> operator*(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    Matrix<T> tmp{lhs};
    tmp *= rhs;
    return tmp;
}

template<typename T = double>
std::vector<T> operator*(const Matrix<T>& lhs, const std::vector<T>& rhs) {
    if(lhs.row_size() != rhs.size())
        throw std::invalid_argument("Matrix*: wrong sizes");

    std::vector<T> tmp(lhs.column_size());
    for(size_t i = 0; i < tmp.size(); ++i)
        for(size_t j = 0; j < lhs.row_size(); ++j)
            tmp[i] += (lhs(i, j) * rhs[j]);

    return tmp;
}

template<typename T = double>
std::vector<T> operator*(const std::vector<T>& lhs, const Matrix<T>& rhs) {
    if(rhs.column_size() != lhs.size())
        throw std::invalid_argument("Matrix*: wrong sizes");

    std::vector<T> tmp(rhs.row_size());
    for(size_t i = 0; i < tmp.size(); ++i)
        for(size_t j = 0; j < lhs.size(); ++j)
            tmp[i] += (lhs[j] * lhs(j, i));

    return tmp;
}

template<typename T = double>
Matrix<T> operator*(T lhs, const Matrix<T>& rhs) {
    Matrix<T> tmp(rhs.column_size(), rhs.row_size());
    for(size_t i = 0; i < tmp.column_size(); ++i)
        for(size_t j = 0; j < tmp.row_size(); ++j)
            tmp(i, j) = lhs * rhs(i, j);

    return tmp;
}

template<typename T = double>
Matrix<T> operator*(const Matrix<T>& lhs, T rhs) {
    return rhs * lhs;
}

template<typename T = double>
class Square_Matrix : public Matrix<T> {
    using Base = Matrix<T>;

public:
    Square_Matrix(size_t size) : Base(size, size) {}
    virtual ~Square_Matrix() = default;

    Square_Matrix(const Lower_Triangle_Matrix<T>& L) : Base(L.size(), L.size()) {
        for(size_t i = 0; i < L.size(); ++i)
            for(size_t j = 0; j <= i; ++j)
                Base(i, j) = L(i, j);
    }

    size_t size() const {
        return Base::column_size();
    }

    static Square_Matrix<T> E(size_t size) {
        Square_Matrix<T> e(size);
        for(size_t i = 0; i < size; ++i)
            e(i, i) = 1.0;
        return e;
    }

    Square_Matrix<T>& operator*=(const Square_Matrix<T>& rhs) {
        if(size() != rhs.size())
            throw std::invalid_argument("Square_Matrix*: wrong sizes");
        Square_Matrix<T> A_new(size());

        for(size_t i = 0; i < A_new.size(); ++i)
            for(size_t j = 0; j < A_new.size(); ++j)
                for(size_t k = 0; k < size(); ++k)
                    A_new(i, j) += (Base::operator()(i, k) * rhs(k, j));

        *this = std::move(A_new);
        return *this;
    }
};

template<typename T = double>
Square_Matrix<T> operator*(const Square_Matrix<T>& lhs, const Square_Matrix<T>& rhs) {
    Square_Matrix<T> tmp{lhs};
    tmp *= rhs;
    return tmp;
}

template<typename T = double>
Square_Matrix<T> operator*(T lhs, const Square_Matrix<T>& rhs) {
    Square_Matrix<T> tmp(rhs.size());
    for(size_t i = 0; i < tmp.column_size(); ++i)
        for(size_t j = 0; j < tmp.row_size(); ++j)
            tmp(i, j) = lhs * rhs(i, j);

    return tmp;
}

template<typename T = double>
Square_Matrix<T> operator*(const Square_Matrix<T>& lhs, T rhs) {
    return rhs * lhs;
}

template<typename T = double>
class Lower_Triangle_Matrix {
    static size_t data_size(size_t size) {
        return (size + 1) * size / 2;
    }

protected:
    size_t size_;
    std::vector<T> data_;

public:
    Lower_Triangle_Matrix(size_t size) : size_(size),
        data_(data_size(size), static_cast<T>(0.0)) {}
    virtual ~Lower_Triangle_Matrix() = default;

    T& operator()(size_t i, size_t j) {
        if(i < j) throw std::invalid_argument("LMatrix: this element upper then diagonal");
        return data_.at(data_size(i) + j);
    }
    const T& operator()(size_t i, size_t j) const {
        if(i < j) throw std::invalid_argument("LMatrix: this element upper then diagonal");
        return data_.at(data_size(i + 1) + j);
    }

    size_t size() const {
        return size_;
    }

    static Lower_Triangle_Matrix<T> E(size_t size) {
        Lower_Triangle_Matrix<T> e(size);
        for(size_t i = 0; i < size; ++i)
            e(i, i) = 1.0;
        return e;
    }
};

template<typename T = double>
Square_Matrix<T> Inverse(Lower_Triangle_Matrix<T> L) {
    Square_Matrix<T> I = Square_Matrix<T>::E(L.size());
    for(size_t i = 0; i < L.size(); ++i) {
        if(fabs(L(i, i)) < 1e-10)
            throw std::invalid_argument("Inverse: degenerate matrix");

        I.mult_row_to_k(i, 1.0 / L(i, i));
        L(i, i) = static_cast<T>(1.0);

        for(size_t j = i + 1; j < L.size(); ++j) {
            T k = -1 * L(j, i);
            L(j, i) = static_cast<T>(0.0);
            I.add_row_with_k(i, j, k);
        }
    }
    return I;
}

template<typename T = double>
std::vector<T> operator+(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    if(lhs.size() != rhs.size())
        throw std::invalid_argument("Vector+: different sizes");

    std::vector<T> tmp(lhs.size());
    for(size_t i = 0; i < tmp.size(); ++i)
        tmp[i] = lhs[i] + rhs[i];

    return tmp;
}

template<typename T = double>
std::vector<T> operator-(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    if(lhs.size() != rhs.size())
        throw std::invalid_argument("Vector+: different sizes");

    std::vector<T> tmp(lhs.size());
    for(size_t i = 0; i < tmp.size(); ++i)
        tmp[i] = lhs[i] - rhs[i];

    return tmp;
}

template<typename T = double>
T N1(const std::vector<T>& v) {
    T n = static_cast<T>(0.0);
    std::for_each(v.begin(), v.end(), [&n](const auto& el){
        if(n < fabs(el)) n = fabs(el);
    });
    return n;
}
