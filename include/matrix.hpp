#pragma once

#include <concepts>

#ifdef QT
#include <QDebug>
#include <QSpan>
#include <QVector>
#else
#include <iostream>
#include <span>
#include <vector>
#endif


template<class T>
concept NumberType = std::integral<T> || std::floating_point<T>;

namespace dynamic_matrix
{
	
class Exception : public std::runtime_error
{
public:
	Exception(const std::string &message) : std::runtime_error(message) {};
};
	
template<NumberType T> using Element = T;
#ifdef QT
template<NumberType T> using Row = QVector<Element<T>>;
template<NumberType T> using Data = QVector<Row<T>>;
template<NumberType T> using InitializerArray = QSpan<const double>;
#else
template<NumberType T> using Row = std::vector<Element<T>>;
template<NumberType T> using Data = std::vector<Row<T>>;
template<NumberType T> using InitializerArray = std::span<const double>;
#endif

template<NumberType T>
class Matrix;
template<NumberType T>
std::ostream& operator<<(std::ostream &os, const Matrix<T> &matrix);
template<NumberType T>
Matrix<T> operator*(const Matrix<T> &mat1, const Matrix<T> &mat2);

template<NumberType T>
class Matrix
{
    template<NumberType K>
    friend std::ostream& operator<<(std::ostream &os, const Matrix<K> &matrix);
    template<NumberType K>
    friend Matrix<K> operator*(const Matrix<K> &mat1, const Matrix<K> &mat2);

public:
    Matrix() noexcept = default;
    Matrix(const std::size_t M, const std::size_t N) noexcept;
    Matrix(const T value, const std::size_t M, const std::size_t N) noexcept;
    Matrix(InitializerArray<T> array, const std::size_t M, const std::size_t N);
	virtual ~Matrix() noexcept = default;

    Matrix(const Matrix<T> &other) noexcept = default;
    Matrix(Matrix<T> &&other) noexcept = default;
    Matrix &operator=(const Matrix<T> &other) noexcept = default;
    Matrix &operator=(Matrix<T> &&other) noexcept = default;

    [[nodiscard]] Matrix operator+(const Matrix<T> &other) noexcept;
    [[nodiscard]] Matrix operator-(const Matrix<T> &other) noexcept;
	[[nodiscard]] Matrix operator*(const T value) noexcept;
	[[nodiscard]] T &operator()(const std::size_t m, const std::size_t n) noexcept;
	[[nodiscard]] const T &operator()(const std::size_t m, const std::size_t n) const noexcept;

	[[nodiscard]] T at(const std::size_t m, const std::size_t n) const;
    [[nodiscard]] Row<T> column(const std::size_t n) const;
    [[nodiscard]] inline const std::size_t &columnCount() const noexcept { return _N_; }
    [[nodiscard]] Data<T> &data();
    [[nodiscard]] const Data<T> &data() const;
	void fill(const T value);
	inline void ones() { fill(static_cast<T>(1)); }
    [[nodiscard]] Row<T> &row(const std::size_t m);
    [[nodiscard]] inline const std::size_t &rowCount() const noexcept { return _M_; }
    [[nodiscard]] const Row<T> &row(const std::size_t m) const;
    void setData(const Data<T> &data) noexcept;
    [[nodiscard]] Matrix<T> transposed() const noexcept;
	inline void zeros() { fill(static_cast<T>(0)); }

protected:
    std::size_t _M_;
    std::size_t _N_;
    Data<T> data_;

};

template<NumberType T>
std::ostream& operator<<(std::ostream &os, const Matrix<T> &matrix)
{
    os << "Matrix(" << matrix._M_ << ", " << matrix._N_ << ")\n{\n";

    for (auto m = 0; m < matrix._M_; ++m)
	{
		os << "\t";
        for (auto n = 0; n < matrix._N_; ++n)
			os << matrix(m, n) << " ";
		os << "\n";
	}

	return os << "}" << std::endl;
}

template<NumberType T>
Matrix<T> operator*(const Matrix<T> &mat1, const Matrix<T> &mat2)
{
    if (mat1._N_ != mat2._M_)
        throw Exception("Column number of first matrix must be equal to row number of second matrix");

    Matrix<T> result(mat1._M_, mat2._N_);

    for (auto m = 0ull; m < result._M_; ++m)
        for (auto n = 0ull; n < result._N_; ++n)
            for (auto r = 0ull; r < mat1._N_; ++r)
				result(m, n) += mat1(m, r) * mat2(r, n);

	return result;
}

template<NumberType T>
Matrix<T>::Matrix(const std::size_t M, const std::size_t N) noexcept
    : _M_(M)
    , _N_(N)
    , data_(_M_, Row<T>(_N_))
{}

template<NumberType T>
Matrix<T>::Matrix(const T value, const std::size_t M, const std::size_t N) noexcept
    : Matrix(M, N)
{
	fill(value);
}

template<NumberType T>
Matrix<T>::Matrix(InitializerArray<T> array, const std::size_t M, const std::size_t N)
    : _M_(M)
    , _N_(N)
{
	if (array.size() < N * M)
		throw Exception("Initialize array is smaller than nessesary");

	for (auto m = 0ull; m < M; ++m)
		for (auto n = 0ull; n < N; ++n)
			data_[m][n] = array[n + m * N];
}

template<NumberType T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) noexcept
{
    Matrix<T> result(other);

	for (auto it1 = result.data_.begin(), it2 = data_.cbegin(); it1 != result.data_.end(); ++it1, ++it2)
		for (auto jt1 = it1->begin(), jt2 = it2->cbegin(); jt1 != it1->end(); ++jt1, ++jt2)
			(*jt1) += (*jt2);

	return result;
}

template<NumberType T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &other) noexcept
{
    Matrix<T> result(other);

	for (auto it1 = result.data_.begin(), it2 = data_.cbegin(); it1 != result.data_.end(); ++it1, ++it2)
		for (auto jt1 = it1->begin(), jt2 = it2->cbegin(); jt1 != it1->end(); ++jt1, ++jt2)
			(*jt1) -= (*it1);

	return result;
}

template<NumberType T>
Matrix<T> Matrix<T>::operator*(const T value) noexcept
{
    Matrix<T> result(*this);

	for (auto it1 = result.data_.begin(); it1 != result.data_.end(); ++it1)
		for (auto jt1 = it1->begin(); jt1 != it1->end(); ++jt1)
			(*jt1) *= value;

	return result;
}

template<NumberType T>
T &Matrix<T>::operator()(const std::size_t m, const std::size_t n) noexcept
{
	return data_[m][n];
}

template<NumberType T>
const T &Matrix<T>::operator()(const std::size_t m, const std::size_t n) const noexcept
{
	return const_cast<const T &>(data_[m][n]);
}

template<NumberType T>
T Matrix<T>::at(const std::size_t m, const std::size_t n) const
{
    if (m >= _M_ || n >= _N_)
		throw Exception("Index out of range");

	return data_[m][n];
}

template<NumberType T>
Row<T> Matrix<T>::column(const std::size_t n) const
{
    if (n >= _N_)
		throw Exception("Index out of range");

    Row<T> result(_M_);

    for (auto m = 0ull; m < _M_; ++m)
		result[m] = data_[m][n];

	return result;
}

template<NumberType T>
Data<T> &Matrix<T>::data()
{
	return data_;
}

template<NumberType T>
const Data<T> &Matrix<T>::data() const
{
    return static_cast<const Data<T> &>(data_);
}

template<NumberType T>
inline void Matrix<T>::fill(const T value)
{
	for (auto &row : data_)
		for (auto &element : row)
			element = value;
}

template<NumberType T>
Row<T> &Matrix<T>::row(const std::size_t m)
{
    if (m >= _M_)
		throw Exception("Index out of range");

	return data_[m];
}

template<NumberType T>
const Row<T> &Matrix<T>::row(const std::size_t m) const
{
    if (m >= _M_)
		throw Exception("Index out of range");

    return static_cast<const Row<T> &>(data_[m]);
}

template<NumberType T>
void Matrix<T>::setData(const Data<T> &data) noexcept
{
    data_ = data;
}

template<NumberType T>
Matrix<T> Matrix<T>::transposed() const noexcept
{
    Matrix<T> result(_M_, _N_);

    for (auto m = 0ull; m < _M_; ++m)
        for (auto n = 0ull; n < _N_; ++n)
			result(n, m) = operator()(m, n);

	return result;
}

}
