#pragma once

#include <Vector.h>

namespace freeaml
{
/**
 * @brief @c Matrix<T> is a class suited for representing dense matrices.
 *
 * This class stores a matrix of elements of type @c T. It overloads the
 * addition (+), subtraction (-), multiplication (*) and division (/) operators
 * for supporting common matrix operations such as matrix addition, matrix
 * multiplication by scalar and matrix multiplication with another matrix or a
 * vector.
 *
 * Some commonly used mathematical operations are also provided in the class
 * (e.g. determining the transpose of the matrix).
 *
 * Support for OpenMP was added to the functions and operators which showed a
 * significant speedup when implemented using multiple threads.
 */
template<typename T>
class Matrix
{
public:
    using size_type = typename Vector<T>::size_type;
    using reference = typename Vector<T>::reference;
    using const_reference = typename Vector<T>::const_reference;

    /** @brief Constructs a matrix with no elements. */
    Matrix();

    /**
     * @brief Constructs a matrix with the contents of an initializer list of
     *        equally-sized initializer lists, each one representing a row of
     *        the matrix.
     * @param init An initializer list holding initializer lists of elements of
     *        type @c T.
     */
    Matrix(std::initializer_list<std::initializer_list<T>> init);

    /**
     * @brief Constructs a matrix with the elements of a vector.
     * @param rows The number of matrix rows.
     * @param cols The number of matrix columns.
     * @param elements A vector holding the <tt>rows × cols</tt> matrix elements
     *        in row-major order (i.e., the elements on the first matrix row
     *        followed by the elements on the second matrix row and so on).
     */
    Matrix(size_type rows, size_type cols, const Vector<T>& elements);

    /**
     * @brief Constructs a matrix with the elements of a vector.
     * @param rows The number of matrix rows.
     * @param cols The number of matrix columns.
     * @param elements A vector holding the <tt>rows × cols</tt> matrix elements
     *        in row-major order (i.e., the elements on the first matrix row
     *        followed by the elements on the second matrix row and so on).
     */
    Matrix(size_type rows, size_type cols, Vector<T>&& elements);

    /**
     * @brief Constructs a matrix with all elements initialized with a value.
     * @param rows The number of matrix rows.
     * @param cols The number of matrix columns.
     * @param x The initializing value for every element of the matrix.
     */
    Matrix(size_type rows, size_type cols, const T& x = T{});

    /**
     * @brief Copy constructor.
     * @param M The matrix from which all elements will be copied.
     */
    Matrix(const Matrix& M) = default;

    /**
     * @brief Move constructor.
     * @param M The matrix from which all elements will be moved.
     */
    Matrix(Matrix&& M) = default;

    /**
     * @brief Returns a reference to a matrix element.
     * @param i The row of the matrix element.
     * @param j The column of the matrix element.
     * @return A reference to the element <tt>(i,j)</tt> of the matrix.
     */
    reference operator()(size_type i, size_type j);

    /**
     * @brief Returns a const reference to a matrix element.
     * @param i The row of the matrix element.
     * @param j The column of the matrix element.
     * @return A const reference to the element <tt>(i,j)</tt> of the matrix.
     */
    const_reference operator()(size_type i, size_type j) const;

    /**
     * @brief Copy-assignment operator.
     * @param M The matrix from which all elements will be copied.
     * @return A reference to @c *this.
     */
    Matrix& operator=(const Matrix& M) = default;

    /**
     * @brief Move-assignment operator.
     * @param M The matrix from which all elements will be moved.
     * @return A reference to @c *this.
     */
    Matrix& operator=(Matrix&& M) = default;

    /**
     * @brief Equality-comparison operator.
     * @param M A matrix to compare against.
     * @return @c true if the matrix is equal to @c M, @c false otherwise.
     */
    bool operator==(const Matrix& M) const;

    /**
     * @brief Inequality-comparison operator.
     * @param M A matrix to compare against.
     * @return @c true if the matrix is not equal to @c M, @c false otherwise.
     */
    bool operator!=(const Matrix& M) const;

    /**
     * @brief Multiplies all elements of the matrix by a scalar.
     * @param c A scalar.
     * @return A reference to @c *this.
     */
    Matrix& operator*=(const T& c);

    /**
     * @brief Divides all elements of the matrix by a scalar.
     * @param c A scalar.
     * @return A reference to @c *this.
     */
    Matrix& operator/=(const T& c);

    /**
     * @brief Performs element-wise addition-assignment with another matrix.
     * @param M A matrix.
     * @return A reference to @c *this.
     */
    Matrix& operator+=(const Matrix& M);

    /**
     * @brief Performs element-wise subtraction-assignment with another matrix.
     * @param M A matrix.
     * @return A reference to @c *this.
     */
    Matrix& operator-=(const Matrix& M);

    /**
     * @brief Computes the transpose of the matrix.
     * @return A copy of the transpose of the matrix.
     */
    Matrix transpose() const;

    /**
     * @brief Checks if the matrix has the same number of rows and columns.
     * @return @c true if the matrix is square, @c false otherwise.
     */
    bool is_square() const;

    /**
     * @brief Checks if the matrix is symmetric.
     * @return @c true if the matrix is symmetric, @c false otherwise.
     */
    bool is_symmetric() const;

    /**
     * @brief Gets the number of rows in the matrix.
     * @return The number of rows in the matrix.
     */
    size_type num_rows() const;

    /**
     * @brief Gets the number of columns in the matrix.
     * @return The number of columns in the matrix.
     */
    size_type num_cols() const;

    /**
     * @brief Checks if the matrix is empty.
     * @return @c true if the matrix is empty, @c false otherwise.
     */
    bool empty() const;

    /**
     * @brief Resizes the matrix.
     * @param rows The new number of matrix rows.
     * @param cols The new number of matrix columns.
     * @param x The initializing value for new elements of the matrix.
     */
    void resize(size_type rows, size_type cols, const T& x = T{});

    /**
     * @brief Clears the matrix.
     */
    void clear();

    /**
     * @brief Returns the matrix elements as a vector (in row-major order).
     * @return The elements of the matrix stored on a vector, with the first
     *         row elements appearing first, then the second row elements and
     *         so on.
     */
    const Vector<T>& flatten() const;

private:
    size_type rows_;     /* number of matrix rows */
    size_type cols_;     /* number of matrix columns */
    Vector<T> elements_; /* matrix elements in row-major order */

}; /* class Matrix<T> */

/**
 * @brief Computes the multiplication of a matrix by a scalar on the right.
 * @param M A matrix.
 * @param c A scalar.
 * @return A copy of @c M with all elements multiplied by @c c.
 */
template<typename T>
Matrix<T> operator*(const Matrix<T>& M, const T& c);

/**
 * @brief Computes the multiplication of a matrix by a scalar on the left.
 * @param c A scalar.
 * @param M A matrix.
 * @return A copy of @c M with all elements multiplied by @c c.
 */
template<typename T>
Matrix<T> operator*(const T& c, const Matrix<T>& M);

/**
 * @brief Computes the multiplication of two matrices.
 * @param M1 a matrix.
 * @param M2 a matrix.
 * @return A matrix which is the result of multipying @c M1 and @c M2.
 */
template<typename T>
Matrix<T> operator*(const Matrix<T>& M1, const Matrix<T>& M2);

/**
 * @brief Computes the multiplication of a matrix and a vector.
 * @param M a matrix.
 * @param v a vector (interpreted as a column vector).
 * @return A vector which is the result of multipying @c M and @c v.
 */
template<typename T>
Vector<T> operator*(const Matrix<T>& M, const Vector<T>& v);

/**
 * @brief Computes the multiplication of a vector and a matrix.
 * @param v a vector (interpreted as a row vector).
 * @param M a matrix.
 * @return A vector which is the result of multipying @c v and @c M.
 */
template<typename T>
Vector<T> operator*(const Vector<T>& v, const Matrix<T>& M);

/**
 * @brief Computes the division of a matrix by a scalar.
 * @param M A matrix.
 * @param c A scalar.
 * @return A copy of @c M with all elements divided by @c c.
 */
template<typename T>
Matrix<T> operator/(const Matrix<T>& M, const T& c);

/**
 * @brief Computes the matrix addition of two equally-sized matrices.
 * @param M1 A matrix.
 * @param M2 A matrix.
 * @return The element-wise sum of @c M1 and @c M2.
 */
template<typename T>
Matrix<T> operator+(const Matrix<T>& M1, const Matrix<T>& M2);

/**
 * @brief Computes the matrix difference of two equally-sized matrices.
 * @param M1 A matrix.
 * @param M2 A matrix.
 * @return The element-wise difference between @c M1 and @c M2.
 */
template<typename T>
Matrix<T> operator-(const Matrix<T>& M1, const Matrix<T>& M2);

/**
 * @brief Computes the element-wise negation of a matrix.
 * @param M A matrix.
 * @return The element-wise negation of @c M.
 */
template<typename T>
Matrix<T> operator-(const Matrix<T>& M);

/**
 * @brief Prints the elements of a matrix to an output stream.
 * @param stream An output stream.
 * @param M A matrix.
 * @return A reference to @c stream.
 */
template<typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T>& M);

/**
 * @brief Generates a random matrix with elements within a given range.
 * @param rows The number of matrix rows.
 * @param cols The number of matrix columns.
 * @param lower_bound The lower bound for the sample interval.
 * @param upper_bound The upper bound for the sample interval.
 * @return A <tt>rows × cols</tt> matrix with elements sampled uniformly from
 *         <tt>[lower_bound, upper_bound]</tt>.
 * @note This function was designed to work only with primitive integer and
 *       floating-point types (e.g. @c int, @c float, @c double etc.).
 */
template<typename T>
Matrix<T> random_matrix(typename Matrix<T>::size_type rows,
                        typename Matrix<T>::size_type cols,
                        const T& lower_bound = T{0},
                        const T& upper_bound = T{1});

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename T>
Matrix<T>::Matrix() : rows_{0}, cols_{0}
{
    /* nothing needs to be done here */
}

template<typename T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> init)
    : rows_(init.size()), cols_(init.size() > 0 ? init.begin()->size() : 0)
{
    if (cols_ == 0)
    {
        clear();
        return;
    }

    elements_.reserve(rows_ * cols_);

    for (const auto& row : init)
    {
        FREEAML_ASSERT(row.size() == cols_);

        for (const T& element : row)
        {
            elements_.push_back(element);
        }
    }
}

template<typename T>
Matrix<T>::Matrix(const size_type rows,
                  const size_type cols,
                  const T& x /* = T{} */)
    : rows_(rows), cols_(cols), elements_(rows * cols, x)
{
    if (rows_ == 0 || cols_ == 0)
    {
        clear();
    }
}

template<typename T>
Matrix<T>::Matrix(const size_type rows,
                  const size_type cols,
                  Vector<T>&& elements)
    : rows_(rows), cols_(cols), elements_(std::move(elements))
{
    FREEAML_ASSERT(rows_ * cols_ == elements_.size());

    if (rows_ == 0 || cols_ == 0)
    {
        clear();
    }
}

template<typename T>
Matrix<T>::Matrix(const size_type rows,
                  const size_type cols,
                  const Vector<T>& elements)
    : rows_(rows), cols_(cols), elements_(elements)
{
    FREEAML_ASSERT(rows_ * cols_ == elements_.size());

    if (rows_ == 0 || cols_ == 0)
    {
        clear();
    }
}

template<typename T>
typename Matrix<T>::reference Matrix<T>::operator()(const size_type i,
                                                    const size_type j)
{
    FREEAML_ASSERT(i < num_rows() && j < num_cols());

    return elements_[i * num_cols() + j];
}

template<typename T>
typename Matrix<T>::const_reference Matrix<T>::operator()(
    const size_type i, const size_type j) const
{
    FREEAML_ASSERT(i < num_rows() && j < num_cols());

    return elements_[i * num_cols() + j];
}

template<typename T>
bool Matrix<T>::operator==(const Matrix<T>& M) const
{
    return rows_ == M.rows_ && cols_ == M.cols_ && elements_ == M.elements_;
}

template<typename T>
bool Matrix<T>::operator!=(const Matrix<T>& M) const
{
    return !operator==(M);
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const T& c)
{
    elements_ *= c;
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator/=(const T& c)
{
    elements_ /= c;
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& M)
{
    FREEAML_ASSERT(num_rows() == M.num_rows());
    FREEAML_ASSERT(num_cols() == M.num_cols());

    elements_ += M.elements_;
    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& M)
{
    FREEAML_ASSERT(num_rows() == M.num_rows());
    FREEAML_ASSERT(num_cols() == M.num_cols());

    elements_ -= M.elements_;
    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::transpose() const
{
    Matrix<T> result(num_cols(), num_rows(), T{});

#ifdef _OPENMP
#pragma omp parallel for
#endif /* _OPENMP */
    for (size_type i = 0; i < num_rows(); ++i)
    {
        for (size_type j = 0; j < num_cols(); ++j)
        {
            result(j, i) = (*this)(i, j);
        }
    }

    return result;
}

template<typename T>
bool Matrix<T>::is_square() const
{
    return num_rows() == num_cols();
}

template<typename T>
bool Matrix<T>::is_symmetric() const
{
    if (is_square() == false)
    {
        return false;
    }

    for (size_type i = 0; i < num_rows(); ++i)
    {
        for (size_type j = i + 1; j < num_cols(); ++j)
        {
            if ((*this)(i, j) != (*this)(j, i))
            {
                return false;
            }
        }
    }

    return true;
}

template<typename T>
typename Matrix<T>::size_type Matrix<T>::num_rows() const
{
    return rows_;
}

template<typename T>
typename Matrix<T>::size_type Matrix<T>::num_cols() const
{
    return cols_;
}

template<typename T>
bool Matrix<T>::empty() const
{
    return elements_.empty();
}

template<typename T>
void Matrix<T>::resize(const size_type rows,
                       const size_type cols,
                       const T& x /* = T{} */)
{
    /* if there is no need for resizing, do nothing */
    if (rows == num_rows() && cols == num_cols())
    {
        return;
    }

    if (rows == 0 || cols == 0)
    {
        clear();
        return;
    }

    Vector<T> elements;
    elements.reserve(rows * cols);

    for (size_type i = 0; i < rows; ++i)
    {
        for (size_type j = 0; j < cols; ++j)
        {
            if (i < num_rows() && j < num_cols())
            {
                elements.push_back((*this)(i, j));
            }
            else
            {
                elements.push_back(x);
            }
        }
    }

    rows_ = rows;
    cols_ = cols;
    elements_ = std::move(elements);
}

template<typename T>
void Matrix<T>::clear()
{
    rows_ = 0;
    cols_ = 0;
    elements_.clear();
}

template<typename T>
const Vector<T>& Matrix<T>::flatten() const
{
    return elements_;
}

template<typename T>
Matrix<T> operator*(const Matrix<T>& M, const T& c)
{
    Matrix<T> result = M;
    result *= c;
    return result;
}

template<typename T>
Matrix<T> operator*(const T& c, const Matrix<T>& M)
{
    return M * c;
}

template<typename T>
Matrix<T> operator*(const Matrix<T>& M1, const Matrix<T>& M2)
{
    FREEAML_ASSERT(M1.num_cols() == M2.num_rows());

    using size_type = typename Matrix<T>::size_type;

    Matrix<T> result(M1.num_rows(), M2.num_cols(), T{});

#ifdef _OPENMP
#pragma omp parallel for
#endif /* _OPENMP */
    for (size_type i = 0; i < result.num_rows(); ++i)
    {
        for (size_type k = 0; k < M1.num_cols(); ++k)
        {
            for (size_type j = 0; j < result.num_cols(); ++j)
            {
                result(i, j) += M1(i, k) * M2(k, j);
            }
        }
    }

    return result;
}

template<typename T>
Vector<T> operator*(const Matrix<T>& M, const Vector<T>& v)
{
    FREEAML_ASSERT(M.num_cols() == v.size());

    using size_type = typename Matrix<T>::size_type;

    Vector<T> result(M.num_rows(), T{});

#ifdef _OPENMP
#pragma omp parallel for
#endif /* _OPENMP */
    for (size_type i = 0; i < result.size(); ++i)
    {
        for (size_type j = 0; j < M.num_cols(); ++j)
        {
            result[i] += M(i, j) * v[j];
        }
    }

    return result;
}

template<typename T>
Vector<T> operator*(const Vector<T>& v, const Matrix<T>& M)
{
    FREEAML_ASSERT(v.size() == M.num_rows());

    using size_type = typename Matrix<T>::size_type;

    Vector<T> result(M.num_cols(), T{});

#ifdef _OPENMP
#pragma omp parallel for
#endif /* _OPENMP */
    for (size_type j = 0; j < result.size(); ++j)
    {
        for (size_type i = 0; i < v.size(); ++i)
        {
            result[j] += v[i] * M(i, j);
        }
    }

    return result;
}

template<typename T>
Matrix<T> operator/(const Matrix<T>& M, const T& c)
{
    Matrix<T> result = M;
    result /= c;
    return result;
}

template<typename T>
Matrix<T> operator+(const Matrix<T>& M1, const Matrix<T>& M2)
{
    Matrix<T> result = M1;
    result += M2;
    return result;
}

template<typename T>
Matrix<T> operator-(const Matrix<T>& M1, const Matrix<T>& M2)
{
    Matrix<T> result = M1;
    result -= M2;
    return result;
}

template<typename T>
Matrix<T> operator-(const Matrix<T>& M)
{
    return Matrix<T>(M.num_rows(), M.num_cols(), -M.flatten());
}

template<typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T>& M)
{
    using size_type = typename Matrix<T>::size_type;

    for (size_type i = 0; i < M.num_rows(); ++i)
    {
        stream << (i == 0 ? "[[" : " [");

        for (size_type j = 0; j < M.num_cols(); ++j)
        {
            stream << M(i, j) << (j + 1 == M.num_cols() ? "" : ", ");
        }

        stream << (i + 1 == M.num_rows() ? "]]\n" : "]\n");
    }

    return stream;
}

template<typename T>
Matrix<T> random_matrix(const typename Matrix<T>::size_type rows,
                        const typename Matrix<T>::size_type cols,
                        const T& lower_bound /* = T{0} */,
                        const T& upper_bound /* = T{1} */)
{
    return Matrix<T>(rows, cols,
                     random_vector<T>(rows * cols, lower_bound, upper_bound));
}

} /* namespace freeaml */
