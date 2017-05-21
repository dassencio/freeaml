#pragma once

#include <QRFactorization.h>

namespace freeaml
{
class QRSolver
{
public:
    /**
     * @brief Solves a linear least-squares problem using QR factorization.
     * @param A The matrix of the linear least-squares problem.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear least-squares problem.
     * @return @c true if a solution could be found, @c false otherwise.
     * @note The number of rows of @c A must be larger than or equal to its
     *       number of columns.
     */
    template<typename MatrixType, typename VectorType>
    static bool solve(const MatrixType& A, VectorType& x, const VectorType& b);

}; /* QRSolver */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename MatrixType, typename VectorType>
bool QRSolver::solve(const MatrixType& A, VectorType& x, const VectorType& b)
{
    using T = typename MatrixType::value_type;
    using size_type = typename MatrixType::size_type;

    /* A is an m × n matrix */
    const size_type n = A.num_cols();

    FREEAML_ASSERT(A.num_rows() >= n);

    Matrix<T> Q, R;

    /* compute the QR factorization of A (A = QR) */
    QRFactorization::factorize(A, Q, R);

    /* let f = (Q^t)*b */
    VectorType f = Q.transpose() * b;

    /*
     * solve R*x = b using backward substitution (R should be upper triangular
     * with nonzero diagonal elements)
     */
    size_type i = n;

    while (i-- > 0)
    {
        T y = T{0};

        for (size_type j = i + 1; j < n; ++j)
        {
            y += R(i, j) * x[j];
        }

        if (R(i, i) == T{0})
        {
            return false;
        }

        x[i] = (f[i] - y) / R(i, i);
    }

    return true;
}

} /* namespace freeaml */
