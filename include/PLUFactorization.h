#pragma once

#include <Matrix.h>
#include <Vector.h>

namespace freeaml
{
class PLUFactorization
{
public:
    /**
     * @brief Constructs a PLU factorizer object.
     * @param enable_partial_pivoting @c true if partial pivoting should be
     *        used, @c false otherwise.
     */
    PLUFactorization(bool enable_partial_pivoting = true);

    /**
     * @brief Computes the PLU factorization of a matrix.
     * @param A The matrix which will be factorized.
     * @param P On success, a permutation matrix.
     * @param L On success, a lower triangular matrix with all diagonal elements
     *        equal to one (1).
     * @param U On success, an nonsingular upper triangular matrix.
     * @return @c true if the factorization <tt>A = PLU</tt> could be built,
     *         @c false otherwise.
     * @note The PLU factorization of @c A exists if and only if @c A is
     *       invertible.
     */
    template<typename MatrixType1,
             typename MatrixType2,
             typename MatrixType3,
             typename MatrixType4>
    bool factorize(const MatrixType1& A,
                   MatrixType2& P,
                   MatrixType3& L,
                   MatrixType4& U) const;

    /**
     * @brief Gets the state set for the partial pivoting technique.
     * @return @c true if partial pivoting is enabled, @c false otherwise.
     */
    bool is_partial_pivoting_enabled() const;

private:
    const bool enable_partial_pivoting_; /* true: enable partial pivoting */

}; /* PLUFactorization */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

PLUFactorization::PLUFactorization(
    const bool enable_partial_pivoting /* = true */)
    : enable_partial_pivoting_(enable_partial_pivoting)
{
    /* nothing needs to be done here */
}

template<typename MatrixType1,
         typename MatrixType2,
         typename MatrixType3,
         typename MatrixType4>
bool PLUFactorization::factorize(const MatrixType1& A,
                                 MatrixType2& P,
                                 MatrixType3& L,
                                 MatrixType4& U) const
{
    FREEAML_ASSERT(A.is_square() == true);

    using T = typename MatrixType1::value_type;
    using size_type = typename MatrixType1::size_type;

    const size_type n = A.num_rows();

    /* use a copy of A to prevent performance issues if A is sparse */
    Matrix<T> K(n, n, A.flatten());

    /* assign the correct dimensions to P, L and U */
    P = MatrixType2(n, n);
    L = MatrixType3(n, n);
    U = MatrixType4(n, n);

    /* perm will be used to build the permutation matrix P at the end */
    Vector<size_type> perm(n);

    /* initialize perm as the "indentity permutation" */
    for (size_type p = 0; p < n; ++p)
    {
        perm[p] = p;
    }

    /* for each row i of K (except for the last one) */
    for (size_type i = 0; i + 1 < n; ++i)
    {
        size_type p = i;

        /*
         * if partial pivoting is enabled:
         *
         *      find p such that |K(p,i)| = max_q|K(perm[q],i)| for i <= q < n
         *
         * if partial pivoting is disabled:
         *
         *      find the first p such that K(perm[p],i) != 0 for i <= p < n
         */
        if (is_partial_pivoting_enabled() == true)
        {
            for (size_type q = i + 1; q < n; ++q)
            {
                if (std::abs(K(perm[p], i)) < std::abs(K(perm[q], i)))
                {
                    p = q;
                }
            }

            /* if K(perm[p],i) = 0, K (A) is not invertible */
            if (K(perm[p], i) == T{0})
            {
                return false;
            }
        }
        else
        {
            while (p < n && K(perm[p], i) == T{0})
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
        std::swap(perm[p], perm[i]);

        for (size_type j = i + 1; j < n; ++j)
        {
            K(perm[j], i) /= K(perm[i], i);

            for (size_type l = i + 1; l < n; ++l)
            {
                K(perm[j], l) -= K(perm[j], i) * K(perm[i], l);
            }
        }
    }

    /* build P */
    for (size_type i = 0; i < n; ++i)
    {
        P(perm[i], i) = T{1};
    }

    /* build L */
    for (size_type i = 0; i < n; ++i)
    {
        L(i, i) = T{1};

        for (size_type j = 0; j < i; ++j)
        {
            L(i, j) = K(perm[i], j);
        }
    }

    /* build U */
    for (size_type i = 0; i < n; ++i)
    {
        for (size_type j = i; j < n; ++j)
        {
            U(i, j) = K(perm[i], j);
        }
    }

    return true;
}

bool PLUFactorization::is_partial_pivoting_enabled() const
{
    return enable_partial_pivoting_;
}

} /* namespace freeaml */
