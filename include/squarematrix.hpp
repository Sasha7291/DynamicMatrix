#pragma once

#include <cmath>

#include "matrix.hpp"


namespace dynamic_matrix
{

template<NumberType T>
class SquareMatrix : public Matrix<T>
{

public:
    SquareMatrix(const std::size_t N) noexcept : Matrix<T>(N, N) {}
    SquareMatrix(const T value, const std::size_t N) noexcept : Matrix<T>(value, N, N) {}
    SquareMatrix(InitializerArray<T> array, const std::size_t N) try : Matrix<T>(array, N, N) {} catch (const Exception &exception) { throw exception; }
    ~SquareMatrix() noexcept = default;

    SquareMatrix(const SquareMatrix<T> &other) noexcept = default;
    SquareMatrix(SquareMatrix<T> &&other) noexcept = default;
    SquareMatrix &operator=(const SquareMatrix<T> &other) noexcept = default;
    SquareMatrix &operator=(SquareMatrix<T> &&other) noexcept = default;

    SquareMatrix(const Matrix<T> &other);
    SquareMatrix(Matrix<T> &&other);
    SquareMatrix &operator=(const Matrix<T> &other);
    SquareMatrix &operator=(Matrix<T> &&other);

    [[nodiscard]] T determinant() const;
    void inverse();
    [[nodiscard]] SquareMatrix inverted() const;
    [[nodiscard]] Matrix<T> toMatrix() const;
    [[nodiscard]] T trace() const;
    [[nodiscard]] SquareMatrix<T> minor(const std::size_t m, const std::size_t n) const;

};

template<NumberType T>
SquareMatrix<T>::SquareMatrix(const Matrix<T> &other)
{
    if (other.rowCount() != other.columnCount())
        throw Exception("Other matrix is not squared");

    this->_M_ = this->_N_ = other.rowCount();
    this->data_ = other.data();
}

template<NumberType T>
SquareMatrix<T>::SquareMatrix(Matrix<T> &&other)
{
    if (other.rowCount() != other.columnCount())
        throw Exception("Other matrix is not squared");

    this->_M_ = this->_N_ = other.rowCount();
    this->data_ = std::move(other.data());
}

template<NumberType T>
SquareMatrix<T> &SquareMatrix<T>::operator=(const Matrix<T> &other)
{
    if (other.rowCount() != other.columnCount())
        throw Exception("Other matrix is not squared");

    this->_M_ = this->_N_ = other.rowCount();
    this->data_ = other.data();
    return *this;
}

template<NumberType T>
SquareMatrix<T> &SquareMatrix<T>::operator=(Matrix<T> &&other)
{
    if (other.rowCount() != other.columnCount())
        throw Exception("Other matrix is not squared");

    this->_M_ = this->_N_ = other.rowCount();
    this->data_ = std::move(other.data());
    return *this;
}

template<NumberType T>
T SquareMatrix<T>::determinant() const
{
    if (this->_N_ == 1ull)
    {
        return this->at(0, 0);
    }
    else if (this->_N_ == 2ull)
    {
        return this->at(0, 0) * this->at(1, 1) - this->at(0, 1) * this->at(1, 0);
    }
    else
    {
        T result = static_cast<T>(0);

        for (auto n = 0ull; n < this->_N_; ++n)
            result += std::pow(-1, n) * this->at(0, n) * minor(0, n).determinant();

        return result;
    }
}

template<NumberType T>
void SquareMatrix<T>::inverse()
{
    SquareMatrix temp(this->transposed());

    double factor = 1.0 / determinant();
    for (auto m = 0ull; m < this->_N_; ++m)
        for (auto n = 0ull; n < this->_N_; ++n)
            this->operator()(m, n) = std::pow(-1, m + n) * factor * temp.minor(m, n).determinant();
}

template<NumberType T>
SquareMatrix<T> SquareMatrix<T>::inverted() const
{
    SquareMatrix result(*this);
    result.inverse();
    return result;
}

template<NumberType T>
Matrix<T> SquareMatrix<T>::toMatrix() const
{
    auto result = Matrix<T>(this->_N_, this->_N_);
    result.setData(this->data_);
    return result;
}

template<NumberType T>
SquareMatrix<T> SquareMatrix<T>::minor(const std::size_t m, const std::size_t n) const
{
    SquareMatrix<T> result(this->_N_ - 1);

    auto row_offset = 0ull;
    auto column_offset = 0ull;
    for (auto i = 0ull; i < this->_N_ - 1; ++i)
    {
        if (i == m)
            ++row_offset;

        for (auto j = 0ull; j < this->_N_ - 1; ++j)
        {
            if (j == n)
                ++column_offset;

            result(i, j) = this->at(i + row_offset, j + column_offset);
        }

        column_offset = 0;
    }

    return result;
}

template<NumberType T>
T SquareMatrix<T>::trace() const
{
    T result = static_cast<T>(0);

    for (auto n = 0; n < this->_N_; ++n)
        result += this->at(n, n);

    return result;
}

}
