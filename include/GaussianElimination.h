#pragma once

#include <IterativeLinearSystemSolverBase.h>
#include <Matrix.h>
#include <Vector.h>
#include <cmath>

namespace freeaml
{
class GaussianElimination
{
public:
    /**
     * @brief Constructs a linear system solver which uses the Gaussian
     *        elimination method.
     * @param enable_partial_pivoting @c true if partial pivoting should be
     *        used, @c false otherwise.
     */
    GaussianElimination(bool enable_partial_pivoting = true);

    /**
     * @brief Solves a linear system using the Gaussian elimination method.
     * @param A The linear system matrix.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear system.
     * @return @c true if the linear system has exactly one solution (i.e., if
     *         @c A is invertible), @c false otherwise.
     */
    template<typename MatrixType, typename VectorType>
    bool solve(const MatrixType& A, VectorType& x, VectorType b);

    /**
     * @brief Gets the state set for the partial pivoting technique.
     * @return @c true if partial pivoting is enabled, @c false otherwise.
     */
    bool is_partial_pivoting_enabled() const;

private:
    const bool enable_partial_pivoting_; /* true: enable partial pivoting */

}; /* class GaussianElimination */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

GaussianElimination::GaussianElimination(
    const bool enable_partial_pivoting /* = true */)
    : enable_partial_pivoting_(enable_partial_pivoting)
{
    /* nothing needs to be done here */
}

template<typename MatrixType, typename VectorType>
bool GaussianElimination::solve(const MatrixType& A,
                                VectorType& x,
                                VectorType b)
{
    using T = typename MatrixType::value_type;
    using size_type = typename MatrixType::size_type;

    FREEAML_ASSERT(
        IterativeLinearSystemSolverBase<T>::check_dimensions(A, x, b) == true);

    const size_type n = A.num_rows();

    /* use a copy of A to prevent performance issues if A is sparse */
    Matrix<T> K(n, n, A.flatten());

    freeaml::Vector<size_type> row(n, 0);

    /* row index vector: used for "swapping" rows of K */
    for (size_type i = 0; i < n; ++i)
    {
        row[i] = i;
    }

    /* for each row i of the matrix K */
    for (size_type i = 0; i < n; ++i)
    {
        size_type p = i;

        /*
         * if partial pivoting is enabled:
         *
         *      find p such that |K(p,i)| = max_q|K(row[q],i)| for i <= q < n
         *
         * if partial pivoting is disabled:
         *
         *      find the first p such that K(row[p],i) != 0 for i <= p < n
         */
        if (is_partial_pivoting_enabled() == true)
        {
            for (size_type q = i + 1; q < n; ++q)
            {
                if (std::abs(K(row[p], i)) < std::abs(K(row[q], i)))
                {
                    p = q;
                }
            }

            /* if K(row[p],i) = 0, K (A) is not invertible */
            if (K(row[p], i) == T{0})
            {
                return false;
            }
        }
        else
        {
            while (p < n && K(row[p], i) == T{0})
            {
                ++p;
            }

            /* if p reaches n, K (A) is not invertible */
            if (p == n)
            {
                return false;
            }
        }

        /* "swap" rows i and p of K */
        std::swap(row[p], row[i]);

        for (size_type j = i + 1; j < n; ++j)
        {
            T m = K(row[j], i) / K(row[i], i);

            for (size_type l = i + 1; l < n; ++l)
            {
                K(row[j], l) -= m * K(row[i], l);
            }

            b[row[j]] -= m * b[row[i]];
        }
    }

    size_type i = n;

    /* backward substitution */
    while (i > 0)
    {
        --i;

        T y = T{0};

        for (size_type j = i + 1; j < n; ++j)
        {
            y += K(row[i], j) * x[j];
        }

        x[i] = (b[row[i]] - y) / K(row[i], i);
    }

    return true;
}

bool GaussianElimination::is_partial_pivoting_enabled() const
{
    return enable_partial_pivoting_;
}

} /* namespace freeaml */
