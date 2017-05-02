#pragma once

#include <Matrix.h>
#include <algorithm>
#include <utility>

namespace freeaml
{
/**
 * @brief @c SparseMatrix<T> is a class suited for representing sparse matrices.
 *
 * This class stores a matrix of elements of type @c T. It overloads the
 * addition (+), subtraction (-), multiplication (*) and division (/) operators
 * for supporting common matrix operations such as matrix addition, matrix
 * multiplication by scalar and matrix multiplication with another matrix or a
 * vector.
 *
 * Whenever possible, this class avoids storing elements equal to @c T{} because
 * @c T{} is zero for all primitive numeric types in C++ (e.g. <tt>int{} ==
 * 0</tt>, <tt>double{} == 0.0</tt> and so on). Therefore, when an element
 * accessed by the user is not explicitly stored in the matrix, it is known to
 * be equal to @c T{}.
 *
 * If the matrix is a @c const object, accessing an element which is not
 * explicitly stored in it incurs no penalties, but if the matrix is not a
 * @c const object, accessing such an element will cause it to actually be
 * introduced in the matrix with value initially set to @c T{}. For that reason,
 * do not access an element of a @c SparseMatrix<T> object @c M unless you are
 * about to set its value to something other than @c T{} or @c M is @c const.
 * Therefore, if @c M is not @c const, only use expressions such as @c M(i,j)
 * if you are setting this element to something not equal to @c T{} (e.g.
 * <tt>M(i,j) = 2.5;</tt>), otherwise you will cause @c M(i,j) to be explicitly
 * stored as @c T{} in the matrix, which is not only wasteful in terms of memory
 * consumption but will also slow down element access considerably.
 *
 * Some commonly used mathematical operations are also provided in the class
 * (e.g. determining the transpose of the matrix).
 *
 * Support for OpenMP was added to the functions and operators which showed a
 * significant speedup when implemented using multiple threads.
 */
template<typename T>
class SparseMatrix
{
public:
    using size_type = typename Vector<T>::size_type;
    using Element = std::pair<size_type, T>;
    using SparseRow = Vector<Element>;

    /** @brief Constructs a sparse matrix with no elements. */
    SparseMatrix();

    /**
     * @brief Constructs a sparse matrix with the contents of an initializer
     *        list of equally-sized initializer lists, each one representing a
     *        row of the matrix.
     * @param init An initializer list holding initializer lists of elements of
     *        type @c T.
     */
    SparseMatrix(std::initializer_list<std::initializer_list<T>> init);

    /**
     * @brief Constructs a sparse matrix with specified dimensions.
     * @param rows The number of matrix rows.
     * @param cols The number of matrix columns.
     */
    SparseMatrix(size_type rows, size_type cols);

    /**
     * @brief Constructs a sparse matrix with the elements of a vector.
     * @param rows The number of matrix rows.
     * @param cols The number of matrix columns.
     * @param elements A vector holding the <tt>rows × cols</tt> matrix elements
     *        in row-major order (i.e., the elements on the first matrix row
     *        followed by the elements on the second matrix row and so on).
     */
    SparseMatrix(size_type rows, size_type cols, const Vector<T>& elements);

    /**
     * @brief Copy constructor.
     * @param M The matrix from which all elements will be copied.
     */
    SparseMatrix(const SparseMatrix& M) = default;

    /**
     * @brief Move constructor.
     * @param M The matrix from which all elements will be moved.
     */
    SparseMatrix(SparseMatrix&& M) = default;

    /**
     * @brief Returns a reference to a matrix element.
     * @param i The row of the matrix element.
     * @param j The column of the matrix element.
     * @return A reference to the element <tt>(i,j)</tt> of the matrix.
     * @note When an element <tt>(i,j)</tt> which is not explicitly stored in
     *       the matrix is accessed, it will be stored and initially set to
     *       @c T{}.
     */
    T& operator()(size_type i, size_type j);

    /**
     * @brief Returns a const reference to a matrix element.
     * @param i The row of the matrix element.
     * @param j The column of the matrix element.
     * @return A const reference to the element <tt>(i,j)</tt> of the matrix.
     * @note When an element <tt>(i,j)</tt> which is not explicitly stored in
     *       the matrix is accessed, a reference to a statically-stored @c T{}
     *       object will be returned, i.e., the element will not be stored in
     *       the matrix due to the fact that it was accessed.
     */
    const T& operator()(size_type i, size_type j) const;

    /**
     * @brief Copy-assignment operator.
     * @param M The matrix from which all elements will be copied.
     * @return A reference to @c *this.
     */
    SparseMatrix& operator=(const SparseMatrix& M) = default;

    /**
     * @brief Move-assignment operator.
     * @param M The matrix from which all elements will be moved.
     * @return A reference to @c *this.
     */
    SparseMatrix& operator=(SparseMatrix&& M) = default;

    /**
     * @brief Equality-comparison operator.
     * @param M A matrix to compare against.
     * @return @c true if the matrix is equal to @c M, @c false otherwise.
     */
    bool operator==(const SparseMatrix& M) const;

    /**
     * @brief Inequality-comparison operator.
     * @param M A matrix to compare against.
     * @return @c true if the matrix is not equal to @c M, @c false otherwise.
     */
    bool operator!=(const SparseMatrix& M) const;

    /**
     * @brief Multiplies all elements of the matrix by a scalar.
     * @param c A scalar.
     * @return A reference to @c *this.
     */
    SparseMatrix& operator*=(const T& c);

    /**
     * @brief Divides all elements of the matrix by a scalar.
     * @param c A scalar.
     * @return A reference to @c *this.
     */
    SparseMatrix& operator/=(const T& c);

    /**
     * @brief Performs element-wise addition-assignment with another matrix.
     * @param M A sparse matrix.
     * @return A reference to @c *this.
     */
    SparseMatrix& operator+=(const SparseMatrix& M);

    /**
     * @brief Performs element-wise subtraction-assignment with another matrix.
     * @param M A sparse matrix.
     * @return A reference to @c *this.
     */
    SparseMatrix& operator-=(const SparseMatrix& M);

    /**
     * @brief Computes the transpose of the matrix.
     * @return A copy of the transpose of the matrix.
     */
    SparseMatrix transpose() const;

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
     */
    void resize(size_type rows, size_type cols);

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
    Vector<T> flatten() const;

    /**
     * @brief Returns the number of elements stored explicitly in the matrix.
     * @return The number of elements stored explicitly in the matrix.
     */
    size_type num_stored() const;

    /**
     * @brief Returns a const reference to the representation of a matrix row.
     * @param i The row number.
     * @return The elements stored on the <tt>i</tt>-th matrix row.
     */
    const SparseRow& row(size_type i) const;

private:
    size_type cols_;             /* number of matrix columns */
    Vector<SparseRow> elements_; /* matrix elements (vector of sparse rows) */

}; /* class SparseMatrix<T> */

/**
 * @brief Computes the multiplication of a sparse matrix by a scalar on the
 *        right.
 * @param M A sparse matrix.
 * @param c A scalar.
 * @return A copy of @c M with all elements multiplied by @c c.
 */
template<typename T>
SparseMatrix<T> operator*(const SparseMatrix<T>& M, const T& c);

/**
 * @brief Computes the multiplication of a sparse matrix by a scalar on the
 *        left.
 * @param c A scalar.
 * @param M A sparse matrix.
 * @return A copy of @c M with all elements multiplied by @c c.
 */
template<typename T>
SparseMatrix<T> operator*(const T& c, const SparseMatrix<T>& M);

/**
 * @brief Computes the multiplication of two sparse matrices.
 * @param M1 a sparse matrix.
 * @param M2 a sparse matrix.
 * @return A sparse matrix which is the result of multipying @c M1 and @c M2.
 */
template<typename T>
SparseMatrix<T> operator*(const SparseMatrix<T>& M1, const SparseMatrix<T>& M2);

/**
 * @brief Computes the multiplication of a sparse matrix and a dense matrix.
 * @param M1 a sparse matrix.
 * @param M2 a dense matrix.
 * @return A dense matrix which is the result of multipying @c M1 and @c M2.
 */
template<typename T>
Matrix<T> operator*(const SparseMatrix<T>& M1, const Matrix<T>& M2);

/**
 * @brief Computes the multiplication of a dense matrix and a sparse matrix.
 * @param M1 a dense matrix.
 * @param M2 a sparse matrix.
 * @return A dense matrix which is the result of multipying @c M1 and @c M2.
 */
template<typename T>
Matrix<T> operator*(const Matrix<T>& M1, const SparseMatrix<T>& M2);

/**
 * @brief Computes the multiplication of a sparse matrix and a vector.
 * @param M a sparse matrix.
 * @param v a vector (interpreted as a column vector).
 * @return A vector which is the result of multipying @c M and @c v.
 */
template<typename T>
Vector<T> operator*(const SparseMatrix<T>& M, const Vector<T>& v);

/**
 * @brief Computes the multiplication of a vector and a sparse matrix.
 * @param v a vector (interpreted as a row vector).
 * @param M a sparse matrix.
 * @return A vector which is the result of multipying @c v and @c M.
 */
template<typename T>
Vector<T> operator*(const Vector<T>& v, const SparseMatrix<T>& M);

/**
 * @brief Computes the division of a sparse matrix by a scalar.
 * @param M A sparse matrix.
 * @param c A scalar.
 * @return A copy of @c M with all elements divided by @c c.
 */
template<typename T>
SparseMatrix<T> operator/(const SparseMatrix<T>& M, const T& c);

/**
 * @brief Computes the matrix addition of two equally-sized sparse matrices.
 * @param M1 A sparse matrix.
 * @param M2 A sparse matrix.
 * @return The element-wise sum of @c M1 and @c M2.
 */
template<typename T>
SparseMatrix<T> operator+(const SparseMatrix<T>& M1, const SparseMatrix<T>& M2);

/**
 * @brief Computes the matrix addition of a dense matrix and an equally-sized
 *        sparse matrix.
 * @param M1 A dense matrix.
 * @param M2 A sparse matrix.
 * @return The element-wise sum of @c M1 and @c M2.
 */
template<typename T>
Matrix<T> operator+(const Matrix<T>& M1, const SparseMatrix<T>& M2);

/**
 * @brief Computes the matrix addition of a sparse matrix and an equally-sized
 *        dense matrix.
 * @param M1 A sparse matrix.
 * @param M2 A dense matrix.
 * @return The element-wise sum of @c M1 and @c M2.
 */
template<typename T>
Matrix<T> operator+(const SparseMatrix<T>& M1, const Matrix<T>& M2);

/**
 * @brief Computes the matrix difference of two equally-sized sparse matrices.
 * @param M1 A sparse matrix.
 * @param M2 A sparse matrix.
 * @return The element-wise difference between @c M1 and @c M2.
 */
template<typename T>
SparseMatrix<T> operator-(const SparseMatrix<T>& M1, const SparseMatrix<T>& M2);

/**
 * @brief Computes the matrix difference of a dense matrix and an equally-sized
 *        sparse matrix.
 * @param M1 A dense matrix.
 * @param M2 A sparse matrix.
 * @return The element-wise difference of @c M1 and @c M2.
 */
template<typename T>
Matrix<T> operator-(const Matrix<T>& M1, const SparseMatrix<T>& M2);

/**
 * @brief Computes the matrix difference of a sparse matrix and an equally-sized
 *        dense matrix.
 * @param M1 A sparse matrix.
 * @param M2 A dense matrix.
 * @return The element-wise difference of @c M1 and @c M2.
 */
template<typename T>
Matrix<T> operator-(const SparseMatrix<T>& M1, const Matrix<T>& M2);

/**
 * @brief Computes the element-wise negation of a sparse matrix.
 * @param M A sparse matrix.
 * @return The element-wise negation of @c M.
 */
template<typename T>
SparseMatrix<T> operator-(const SparseMatrix<T>& M);

/**
 * @brief Prints the elements of a sparse matrix to an output stream.
 * @param stream An output stream.
 * @param M A sparse matrix.
 * @return A reference to @c stream.
 */
template<typename T>
std::ostream& operator<<(std::ostream& stream, const SparseMatrix<T>& M);

/**
 * @brief Generates a random sparse matrix with elements within a given range.
 * @param rows The number of matrix rows.
 * @param cols The number of matrix columns.
 * @param nonzero The number of nonzero elements.
 * @param lower_bound The lower bound for the sample interval.
 * @param upper_bound The upper bound for the sample interval.
 * @return A <tt>rows × cols</tt> matrix with @c nonzero elements sampled
 *         uniformly from <tt>[lower_bound, upper_bound]</tt>.
 * @note This function was designed to work only with primitive integer and
 *       floating-point types (e.g. @c int, @c float, @c double etc.).
 */
template<typename T>
SparseMatrix<T> random_sparse_matrix(
    const typename SparseMatrix<T>::size_type rows,
    const typename SparseMatrix<T>::size_type cols,
    const typename SparseMatrix<T>::size_type nonzero,
    const T& lower_bound = T{0},
    const T& upper_bound = T{1});

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename T>
SparseMatrix<T>::SparseMatrix() : cols_{0}, elements_{}
{
    /* nothing needs to be done here */
}

template<typename T>
SparseMatrix<T>::SparseMatrix(
    std::initializer_list<std::initializer_list<T>> init)
    : cols_(init.size() > 0 ? init.begin()->size() : 0), elements_(init.size())
{
    if (cols_ == 0)
    {
        clear();
        return;
    }

    size_type i = 0;

    for (const auto& row : init)
    {
        FREEAML_ASSERT(row.size() == cols_);

        size_type j = 0;

        for (const T& element : row)
        {
            if (element != T{})
            {
                elements_[i].push_back(Element(j, element));
            }
            ++j;
        }
        ++i;
    }
}

template<typename T>
SparseMatrix<T>::SparseMatrix(const size_type rows, const size_type cols)
    : cols_(cols), elements_(rows)
{
    if (rows == 0 || cols == 0)
    {
        clear();
    }
}

template<typename T>
SparseMatrix<T>::SparseMatrix(const size_type rows,
                              const size_type cols,
                              const Vector<T>& elements)
    : cols_(cols), elements_(rows)
{
    if (rows == 0 || cols == 0)
    {
        clear();
        return;
    }

    FREEAML_ASSERT(rows * cols == elements.size());

    for (size_type i = 0; i < rows; ++i)
    {
        for (size_type j = 0; j < cols; ++j)
        {
            const T& element = elements[i * cols + j];

            if (element != T{})
            {
                elements_[i].push_back(Element(j, element));
            }
        }
    }
}

template<typename T>
T& SparseMatrix<T>::operator()(const size_type i, const size_type j)
{
    FREEAML_ASSERT(i < num_rows() && j < num_cols());

    for (Element& element : elements_[i])
    {
        if (element.first == j)
        {
            return element.second;
        }
    }

    elements_[i].push_back(Element(j, T{}));

    return elements_[i].back().second;
}

template<typename T>
const T& SparseMatrix<T>::operator()(const size_type i, const size_type j) const
{
    FREEAML_ASSERT(i < num_rows() && j < num_cols());

    static T zero = T{};

    for (const Element& element : elements_[i])
    {
        if (element.first == j)
        {
            return element.second;
        }
    }

    return zero;
}

template<typename T>
bool SparseMatrix<T>::operator==(const SparseMatrix<T>& M) const
{
    if (num_rows() != M.num_rows() || num_cols() != M.num_cols())
    {
        return false;
    }

    for (size_type i = 0; i < num_rows(); ++i)
    {
        for (const Element& element : elements_[i])
        {
            if (element.second != M(i, element.first))
            {
                return false;
            }
        }
    }

    return true;
}

template<typename T>
bool SparseMatrix<T>::operator!=(const SparseMatrix<T>& M) const
{
    return !operator==(M);
}

template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator*=(const T& c)
{
    for (SparseRow& row : elements_)
    {
        for (Element& element : row)
        {
            element.second *= c;
        }
    }

    return *this;
}

template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator/=(const T& c)
{
    for (SparseRow& row : elements_)
    {
        for (Element& element : row)
        {
            element.second /= c;
        }
    }

    return *this;
}

template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator+=(const SparseMatrix<T>& M)
{
    FREEAML_ASSERT(num_rows() == M.num_rows());
    FREEAML_ASSERT(num_cols() == M.num_cols());

    for (size_type i = 0; i < num_rows(); ++i)
    {
        for (const Element& element : M.elements_[i])
        {
            (*this)(i, element.first) += element.second;
        }
    }

    return *this;
}

template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator-=(const SparseMatrix<T>& M)
{
    FREEAML_ASSERT(num_rows() == M.num_rows());
    FREEAML_ASSERT(num_cols() == M.num_cols());

    for (size_type i = 0; i < num_rows(); ++i)
    {
        for (const Element& element : M.elements_[i])
        {
            (*this)(i, element.first) -= element.second;
        }
    }

    return *this;
}

template<typename T>
SparseMatrix<T> SparseMatrix<T>::transpose() const
{
    SparseMatrix<T> result(num_cols(), num_rows());

    /* operator() on SparseMatrix is not thread safe, so no parallelism here */
    for (size_type i = 0; i < num_rows(); ++i)
    {
        for (const Element& element : elements_[i])
        {
            result(element.first, i) = element.second;
        }
    }

    return result;
}

template<typename T>
bool SparseMatrix<T>::is_square() const
{
    return num_rows() == num_cols();
}

template<typename T>
bool SparseMatrix<T>::is_symmetric() const
{
    if (is_square() == false)
    {
        return false;
    }

    for (size_type i = 0; i < num_rows(); ++i)
    {
        for (const Element& element : elements_[i])
        {
            if ((*this)(element.first, i) != element.second)
            {
                return false;
            }
        }
    }

    return true;
}

template<typename T>
typename SparseMatrix<T>::size_type SparseMatrix<T>::num_rows() const
{
    return elements_.size();
}

template<typename T>
typename SparseMatrix<T>::size_type SparseMatrix<T>::num_cols() const
{
    return cols_;
}

template<typename T>
bool SparseMatrix<T>::empty() const
{
    return elements_.empty();
}

template<typename T>
void SparseMatrix<T>::resize(const size_type rows, const size_type cols)
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

    /* remove elements with invalid column indices (after resizing) */
    if (cols < cols_)
    {
        auto predicate = [cols](const Element& element) {
            return element.first >= cols;
        };

        for (SparseRow& row : elements_)
        {
            row.erase(std::remove_if(row.begin(), row.end(), predicate),
                      row.end());
        }
    }

    cols_ = cols;
    elements_.resize(rows);
}

template<typename T>
void SparseMatrix<T>::clear()
{
    cols_ = 0;
    elements_.clear();
}

template<typename T>
Vector<T> SparseMatrix<T>::flatten() const
{
    Vector<T> elements(num_rows() * num_cols(), T{});

    for (size_type i = 0; i < num_rows(); ++i)
    {
        for (const Element& element : elements_[i])
        {
            elements[i * num_cols() + element.first] = element.second;
        }
    }

    return elements;
}

template<typename T>
typename SparseMatrix<T>::size_type SparseMatrix<T>::num_stored() const
{
    size_type count = 0;

    for (const SparseRow& row : elements_)
    {
        count += row.size();
    }

    return count;
}

template<typename T>
const typename SparseMatrix<T>::SparseRow& SparseMatrix<T>::row(
    const size_type i) const
{
    FREEAML_ASSERT(i < num_rows());

    return elements_[i];
}

template<typename T>
SparseMatrix<T> operator*(const SparseMatrix<T>& M, const T& c)
{
    SparseMatrix<T> result = M;
    result *= c;
    return result;
}

template<typename T>
SparseMatrix<T> operator*(const T& c, const SparseMatrix<T>& M)
{
    return M * c;
}

template<typename T>
SparseMatrix<T> operator*(const SparseMatrix<T>& M1, const SparseMatrix<T>& M2)
{
    FREEAML_ASSERT(M1.num_cols() == M2.num_rows());

    using size_type = typename SparseMatrix<T>::size_type;
    using Element = typename SparseMatrix<T>::Element;

    SparseMatrix<T> result(M1.num_rows(), M2.num_cols());

#ifdef _OPENMP
#pragma omp parallel for
#endif /* _OPENMP */
    for (size_type i = 0; i < result.num_rows(); ++i)
    {
        for (const Element& element1 : M1.row(i))
        {
            for (const Element& element2 : M2.row(element1.first))
            {
                result(i, element2.first) += element1.second * element2.second;
            }
        }
    }

    return result;
}

template<typename T>
Matrix<T> operator*(const SparseMatrix<T>& M1, const Matrix<T>& M2)
{
    FREEAML_ASSERT(M1.num_cols() == M2.num_rows());

    using size_type = typename SparseMatrix<T>::size_type;
    using Element = typename SparseMatrix<T>::Element;

    Matrix<T> result(M1.num_rows(), M2.num_cols(), T{});

#ifdef _OPENMP
#pragma omp parallel for
#endif /* _OPENMP */
    for (size_type i = 0; i < result.num_rows(); ++i)
    {
        for (const Element& element : M1.row(i))
        {
            for (size_type j = 0; j < result.num_cols(); ++j)
            {
                result(i, j) += element.second * M2(element.first, j);
            }
        }
    }

    return result;
}

template<typename T>
Matrix<T> operator*(const Matrix<T>& M1, const SparseMatrix<T>& M2)
{
    FREEAML_ASSERT(M1.num_cols() == M2.num_rows());

    using size_type = typename SparseMatrix<T>::size_type;
    using Element = typename SparseMatrix<T>::Element;

    Matrix<T> result(M1.num_rows(), M2.num_cols(), T{});

#ifdef _OPENMP
#pragma omp parallel for
#endif /* _OPENMP */
    for (size_type i = 0; i < result.num_rows(); ++i)
    {
        for (size_type k = 0; k < M1.num_cols(); ++k)
        {
            for (const Element& element : M2.row(k))
            {
                result(i, element.first) += M1(i, k) * element.second;
            }
        }
    }

    return result;
}

template<typename T>
Vector<T> operator*(const SparseMatrix<T>& M, const Vector<T>& v)
{
    FREEAML_ASSERT(M.num_cols() == v.size());

    using size_type = typename SparseMatrix<T>::size_type;
    using Element = typename SparseMatrix<T>::Element;

    Vector<T> result(M.num_rows(), T{});

#ifdef _OPENMP
#pragma omp parallel for
#endif /* _OPENMP */
    for (size_type i = 0; i < result.size(); ++i)
    {
        for (const Element& element : M.row(i))
        {
            result[i] += element.second * v[element.first];
        }
    }

    return result;
}

template<typename T>
Vector<T> operator*(const Vector<T>& v, const SparseMatrix<T>& M)
{
    FREEAML_ASSERT(v.size() == M.num_rows());

    using size_type = typename SparseMatrix<T>::size_type;
    using Element = typename SparseMatrix<T>::Element;

    Vector<T> result(M.num_cols(), T{});

#ifdef _OPENMP
#pragma omp parallel
    {
        Vector<T> local_result(M.num_cols(), T{});

#pragma omp for nowait
        for (size_type i = 0; i < M.num_rows(); ++i)
        {
            for (const Element& element : M.row(i))
            {
                local_result[element.first] += v[i] * element.second;
            }
        }
#pragma omp critical
        {
            result += local_result;
        }
    }
#else
    /* serial implementation */
    for (size_type i = 0; i < M.num_rows(); ++i)
    {
        for (const Element& element : M.row(i))
        {
            result[element.first] += v[i] * element.second;
        }
    }
#endif /* _OPENMP */

    return result;
}

template<typename T>
SparseMatrix<T> operator/(const SparseMatrix<T>& M, const T& c)
{
    SparseMatrix<T> result = M;
    result /= c;
    return result;
}

template<typename T>
SparseMatrix<T> operator+(const SparseMatrix<T>& M1, const SparseMatrix<T>& M2)
{
    SparseMatrix<T> result = M1;
    result += M2;
    return result;
}

template<typename T>
Matrix<T> operator+(const Matrix<T>& M1, const SparseMatrix<T>& M2)
{
    FREEAML_ASSERT(M1.num_rows() == M2.num_rows());
    FREEAML_ASSERT(M1.num_cols() == M2.num_cols());

    using size_type = typename Matrix<T>::size_type;
    using Element = typename SparseMatrix<T>::Element;

    Matrix<T> result = M1;

    for (size_type i = 0; i < result.num_rows(); ++i)
    {
        for (const Element& element : M2.row(i))
        {
            result(i, element.first) += element.second;
        }
    }

    return result;
}

template<typename T>
Matrix<T> operator+(const SparseMatrix<T>& M1, const Matrix<T>& M2)
{
    return M2 + M1;
}

template<typename T>
SparseMatrix<T> operator-(const SparseMatrix<T>& M1, const SparseMatrix<T>& M2)
{
    SparseMatrix<T> result = M1;
    result -= M2;
    return result;
}

template<typename T>
Matrix<T> operator-(const Matrix<T>& M1, const SparseMatrix<T>& M2)
{
    FREEAML_ASSERT(M1.num_rows() == M2.num_rows());
    FREEAML_ASSERT(M1.num_cols() == M2.num_cols());

    using size_type = typename Matrix<T>::size_type;
    using Element = typename SparseMatrix<T>::Element;

    Matrix<T> result = M1;

    for (size_type i = 0; i < result.num_rows(); ++i)
    {
        for (const Element& element : M2.row(i))
        {
            result(i, element.first) -= element.second;
        }
    }

    return result;
}

template<typename T>
Matrix<T> operator-(const SparseMatrix<T>& M1, const Matrix<T>& M2)
{
    FREEAML_ASSERT(M1.num_rows() == M2.num_rows());
    FREEAML_ASSERT(M1.num_cols() == M2.num_cols());

    using size_type = typename Matrix<T>::size_type;
    using Element = typename SparseMatrix<T>::Element;

    Matrix<T> result = -M2;

    for (size_type i = 0; i < result.num_rows(); ++i)
    {
        for (const Element& element : M1.row(i))
        {
            result(i, element.first) += element.second;
        }
    }

    return result;
}

template<typename T>
SparseMatrix<T> operator-(const SparseMatrix<T>& M)
{
    SparseMatrix<T> result(M.num_rows(), M.num_cols());

    using size_type = typename SparseMatrix<T>::size_type;
    using Element = typename SparseMatrix<T>::Element;

    for (size_type i = 0; i < M.num_rows(); ++i)
    {
        for (const Element& element : M.row(i))
        {
            result(i, element.first) = -element.second;
        }
    }

    return result;
}

template<typename T>
std::ostream& operator<<(std::ostream& stream, const SparseMatrix<T>& M)
{
    using size_type = typename SparseMatrix<T>::size_type;

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
SparseMatrix<T> random_sparse_matrix(
    const typename SparseMatrix<T>::size_type rows,
    const typename SparseMatrix<T>::size_type cols,
    const typename SparseMatrix<T>::size_type nonzero,
    const T& lower_bound /* = T{0} */,
    const T& upper_bound /* = T{1} */)
{
    FREEAML_ASSERT(nonzero <= rows * cols);
    FREEAML_ASSERT(lower_bound < upper_bound);

    if (rows == 0 || cols == 0)
    {
        return {};
    }

    using DistributionType =
        typename std::conditional<std::is_integral<T>::value,
                                  std::uniform_int_distribution<T>,
                                  std::uniform_real_distribution<T>>::type;

    using size_type = typename SparseMatrix<T>::size_type;

    std::random_device device;
    std::mt19937_64 generator(device());

    std::uniform_int_distribution<size_type> row_chooser(0, rows - 1);
    std::uniform_int_distribution<size_type> col_chooser(0, cols - 1);
    DistributionType distribution(lower_bound, upper_bound);

    SparseMatrix<T> result(rows, cols);

    size_type count = 0;

    while (count < nonzero)
    {
        /* choose an element (i,j) */
        size_type i = row_chooser(generator);
        size_type j = col_chooser(generator);

        T x = distribution(generator);

        if (x != T{} && result(i, j) == T{})
        {
            result(i, j) = x;
            ++count;
        }
    }

    return result;
}

} /* namespace freeaml */
