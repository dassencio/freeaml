#pragma once

#include <Matrix.h>
#include <Vector.h>
#include <rotation.h>

namespace freeaml
{
class BidiagonalFactorization
{
public:
    /**
     * @brief Computes the upper bidiagonal factorization of a matrix.
     * @param A The matrix which will be factorized.
     * @param U On success, an orthogonal matrix.
     * @param Z On success, an upper bidiagonal matrix.
     * @param V On success, an orthogonal matrix.
     * @return @c true if the factorization <tt>A = UZV<sup>t</sup></tt> could
     *         be built, @c false otherwise.
     * @note Householder transformations are used to compute @c U, @c V and
     *       @c Z.
     */
    template<typename MatrixType1,
             typename MatrixType2,
             typename MatrixType3,
             typename MatrixType4>
    static bool factorize(const MatrixType1& A,
                          MatrixType2& U,
                          MatrixType3& Z,
                          MatrixType4& V);

}; /* BidiagonalFactorization */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename MatrixType1,
         typename MatrixType2,
         typename MatrixType3,
         typename MatrixType4>
bool BidiagonalFactorization::factorize(const MatrixType1& A,
                                        MatrixType2& U,
                                        MatrixType3& Z,
                                        MatrixType4& V)
{
    using T = typename MatrixType1::value_type;
    using size_type = typename MatrixType1::size_type;

    const size_type m = A.num_rows();
    const size_type n = A.num_cols();

    /* use a copy of A to prevent performance issues if A is sparse */
    Matrix<T> K(m, n, A.flatten());

    /*
     * assign the correct dimensions to U and V (U and V will be products of
     * Householder transformations)
     */
    U = MatrixType2(m, m);
    V = MatrixType4(n, n);

    /* initialize U to the m × m identity matrix */
    for (size_type i = 0; i < m; ++i)
    {
        U(i, i) = T{1};
    }

    /* initialize V to the n × n identity matrix */
    for (size_type i = 0; i < n; ++i)
    {
        V(i, i) = T{1};
    }

    size_type d = std::min(m, n);

    /*
     * for each column p of B which must have the portion below its subdiagonal
     * zeroed out
     */
    for (size_type p = 0; p < d; ++p)
    {
        size_type M = m - p;
        size_type N = n - p;

        /* let uL := K(p:m-1, p) */
        Vector<T> uL(M);

        for (size_type i = 0; i < M; ++i)
        {
            uL[i] = K(i + p, p);
        }

        /*
         * define the Householder transformation for uL which projects it onto
         * e_0 = (1, 0, ..., 0)
         */
        Vector<T> vL;

        T b = T{0};
        T c = T{0};

        householder(uL, vL, b, c, 0);

        /* If b = 0, no Householder transformation is necessary */
        if (b != T{0})
        {
            Matrix<T> B(M, N);
            Matrix<T> C(M, m);

            for (size_type i = 0; i < M; ++i)
            {
                /* B <-- (I - bvv^t)K(p:m-1, p:n-1), where v = vL */
                for (size_type j = 0; j < N; ++j)
                {
                    B(i, j) = K(i + p, j + p);

                    for (size_type k = 0; k < M; ++k)
                    {
                        B(i, j) -= b * vL[i] * vL[k] * K(k + p, j + p);
                    }
                }

                /* C <-- (I - bvv^t)U(p:m-1, 0:m-1), where v = vL */
                for (size_type j = 0; j < m; ++j)
                {
                    C(i, j) = U(i + p, j);

                    for (size_type k = 0; k < M; ++k)
                    {
                        C(i, j) -= b * vL[i] * vL[k] * U(k + p, j);
                    }
                }
            }

            for (size_type i = 0; i < M; ++i)
            {
                /* K(p:m-1, p:n-1) <-- -c*B */
                for (size_type j = 0; j < N; ++j)
                {
                    K(i + p, j + p) = -c * B(i, j);
                }

                /* U(p:m-1, 0:m-1) <-- -c*C */
                for (size_type j = 0; j < m; ++j)
                {
                    U(i + p, j) = -c * C(i, j);
                }
            }
        }

        /* if V needs to be updated (the "last" V is just the identity) */
        if (N > 1)
        {
            /* let uR := K(p, p+1:n-1) */
            Vector<T> uR(N - 1);

            for (size_type j = 0; j < N - 1; ++j)
            {
                uR[j] = K(p, p + j + 1);
            }

            /*
             * define the Householder transformation for uR which projects it
             * onto e_0 = (1, 0, ..., 0)
             */
            Vector<T> vR;

            c = T{0};
            b = T{0};

            householder(uR, vR, b, c, 0);

            /* If b = 0, no Householder transformation is necessary */
            if (b != T{0})
            {
                Matrix<T> D(M, N - 1);
                Matrix<T> E(n, N - 1);

                for (size_type j = 0; j < N - 1; ++j)
                {
                    /* D <-- K(p:m-1, p+1:n-1)(I - b(vv^t)^t), where v = vR */
                    for (size_type i = 0; i < M; ++i)
                    {
                        D(i, j) = K(p + i, p + j + 1);

                        for (size_type k = 0; k < N - 1; ++k)
                        {
                            D(i, j) -= K(p + i, p + k + 1) * b * vR[k] * vR[j];
                        }
                    }

                    /* E = V(0:n-1, p+1:n-1)(I - b(vv^t)^t), where v = vR */
                    for (size_type i = 0; i < n; ++i)
                    {
                        E(i, j) = V(i, p + j + 1);

                        for (size_type k = 0; k < N - 1; ++k)
                        {
                            E(i, j) -= V(i, p + k + 1) * b * vR[k] * vR[j];
                        }
                    }
                }

                for (size_type j = 0; j < N - 1; ++j)
                {
                    /* K(p:m-1, p+1:n-1) <-- -c*D */
                    for (size_type i = 0; i < M; ++i)
                    {
                        K(p + i, p + j + 1) = -c * D(i, j);
                    }

                    /* V(0:n-1, p+1:n-1) <-- -c*E */
                    for (size_type i = 0; i < n; ++i)
                    {
                        V(i, p + j + 1) = -c * E(i, j);
                    }
                }
            }
        }
    }

    /*
     * Z <-- K (at this point, K is the upper biadiagonal form of the original
     * matrix A)
     */
    Z = MatrixType3(m, n);

    for (size_type i = 0; i < d; ++i)
    {
        Z(i, i) = K(i, i);

        if (i + 1 < n)
        {
            Z(i, i + 1) = K(i, i + 1);
        }
    }

    /* U <-- U^t */
    U = U.transpose();

    return true;
}

} /* namespace freeaml */
