#pragma once

#include <cmath>

namespace freeaml
{
class CholeskyFactorization
{
public:
    /**
     * @brief Computes the Cholesky factorization of a positive-definite
     *        symmetric matrix.
     * @param A The matrix which will be factorized.
     * @param L On success, a lower triangular matrix with positive diagonal
     *        entries such that <tt>A = LL<sup>t</sup></tt>.
     * @return @c true if the factorization <tt>A = LL<sup>t</sup></tt> could be
     *         built, @c false otherwise.
     */
    template<typename MatrixType1, typename MatrixType2>
    static bool factorize(const MatrixType1& A, MatrixType2& L);

}; /* CholeskyFactorization */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename MatrixType1, typename MatrixType2>
bool CholeskyFactorization::factorize(const MatrixType1& A, MatrixType2& L)
{
    using T = typename MatrixType1::value_type;
    using size_type = typename MatrixType1::size_type;

    const size_type n = A.num_rows();

    /* assign the correct dimensions to L */
    L = MatrixType2(n, n);

    for (size_type j = 0; j < n; ++j)
    {
        L(j, j) = A(j, j);

        for (size_type k = 0; k < j; ++k)
        {
            L(j, j) -= L(j, k) * L(j, k);
        }

        /*
         * if L(j,j) <= 0, A is not positive definite and therefore has no
         * Cholesky factorization
         */
        if (L(j, j) > T{0})
        {
            L(j, j) = std::sqrt(L(j, j));
        }
        else
        {
            return false;
        }

        for (size_type i = j + 1; i < n; ++i)
        {
            L(i, j) = A(i, j);

            for (size_type k = 0; k < j; ++k)
            {
                L(i, j) -= L(i, k) * L(j, k);
            }

            L(i, j) /= L(j, j);
        }
    }

    return true;
}

} /* namespace freeaml */
