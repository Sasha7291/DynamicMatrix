#pragma once

#include "matrix.hpp"

#include <complex>


namespace dynamic_matrix
{

    template<NumberType T>
    class ComplexMatrix : public Matrix<std::complex<T>>
    {

    public:
        ComplexMatrix(const std::size_t M, const std::size_t N) noexcept : Matrix<T>(M, N) {}
        ComplexMatrix(const std::complex<T> value, const std::size_t M, const std::size_t N) noexcept : Matrix<std::complex<T>>(value, M, N) {}
        ComplexMatrix(std::span<std::complex<T>> array, const std::size_t M, const std::size_t N)
        try : Matrix<std::complex<T>>(array, M, N) {}
        catch (const Exception &exception) { throw exception; }
        ~ComplexMatrix() noexcept = default;

        void conjugate();
        [[nodiscard]] ComplexMatrix conjugated() const;
        void conjugatedTranspose();
        [[nodiscard]] ComplexMatrix conjugatedTransposed() const;

    };

    template<NumberType T>
    void ComplexMatrix<T>::conjugate()
    {
        for (auto m = 0ull; m < this->_M_; ++m)
            for (auto n = 0ull; n < this->_N_; ++n)
                this->operator()(m, n).conj();
    }

    template<NumberType T>
    ComplexMatrix<T> ComplexMatrix<T>::conjugated() const
    {
        ComplexMatrix result(*this);
        result.conjugate();
        return result;
    }

    template<NumberType T>
    void ComplexMatrix<T>::conjugatedTranspose()
    {
        *this = this->transposed();
        conjugate();
    }

    template<NumberType T>
    ComplexMatrix<T> ComplexMatrix<T>::conjugatedTransposed() const
    {
        ComplexMatrix result(*this);
        result.conjugatedTranspose();
        return result;
    }

}
