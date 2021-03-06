#pragma once

#include <Matrix.h>
#include <rotation.h>
#include <cmath>

namespace freeaml
{
class QRFactorization
{
public:
    /**
     * @brief Computes the QR factorization of a matrix.
     * @param A The matrix which will be factorized.
     * @param Q At the end, an orthogonal matrix.
     * @param R At the end, an upper triangular matrix with nonnegative diagonal
     *        elements.
     * @note The number of rows of @c A must be larger than or equal to its
     *       number of columns.
     */
    template<typename MatrixType1, typename MatrixType2, typename MatrixType3>
    static void factorize(const MatrixType1& A, MatrixType2& Q, MatrixType3& R);

}; /* QRFactorization */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename MatrixType1, typename MatrixType2, typename MatrixType3>
void QRFactorization::factorize(const MatrixType1& A,
                                MatrixType2& Q,
                                MatrixType3& R)
{
    using T = typename MatrixType1::value_type;
    using size_type = typename MatrixType1::size_type;

    const size_type m = A.num_rows();
    const size_type n = A.num_cols();

    /* use a copy of A to prevent performance issues if A is sparse */
    Matrix<T> K(m, n, A.flatten());

    /* at the end, P will be the product of all Givens rotations applied */
    Matrix<T> P = identity_matrix<T>(m);

    size_type d = std::min(m, n);

    /* assign the correct dimensions to Q and R */
    Q = MatrixType2(m, d);
    R = MatrixType3(d, n);

    /*
     * for each column p of K which must have its subdiagonal elements zeroed
     * out
     */
    for (size_type p = 0; p < d; ++p)
    {
        /* for each element on the p-th column of K which must be zeroed out */
        for (size_type i = p + 1; i < m; ++i)
        {
            /* if K(i,p) is not already zero */
            if (K(i, p) != T{0})
            {
                /*
                 * define the Givens rotation matrix parameters c and s for
                 * zeroing out K(i,p)
                 */
                T c = T{0};
                T s = T{0};

                givens(K(p, p), K(i, p), c, s);

                /* apply the Givens rotation to K */
                for (size_type j = p; j < n; ++j)
                {
                    T a = K(p, j);
                    T b = K(i, j);

                    K(p, j) = c * a - s * b;
                    K(i, j) = s * a + c * b;
                }

                /* apply the Givens rotation to P */
                for (size_type j = 0; j < m; ++j)
                {
                    T a = P(p, j);
                    T b = P(i, j);

                    P(p, j) = c * a - s * b;
                    P(i, j) = s * a + c * b;
                }
            }
        }
    }

    /* build Q from the first m rows of P^t */
    for (size_type i = 0; i < m; ++i)
    {
        for (size_type j = 0; j < d; ++j)
        {
            Q(i, j) = P(j, i);
        }
    }

    /* at this point, K is upper triangular; build R from its first d rows */
    for (size_type i = 0; i < d; ++i)
    {
        for (size_type j = i; j < n; ++j)
        {
            R(i, j) = K(i, j);
        }
    }
}

} /* namespace freeaml */
