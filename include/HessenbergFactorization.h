#pragma once

#include <Matrix.h>
#include <Vector.h>
#include <rotation.h>

namespace freeaml
{
class HessenbergFactorization
{
public:
    /**
     * @brief Computes the Hessenberg factorization of a matrix.
     * @param A The matrix which will be factorized.
     * @param Q At the end, an orthogonal matrix.
     * @param H At the end, an upper Hessenberg matrix with nonnegative
     *        elements on its first subdiagonal.
     * @note The Hessenberg factorization of @c A is
     *       <tt>A = QHQ<sup>t</sup></tt>.
     */
    template<typename MatrixType1, typename MatrixType2, typename MatrixType3>
    static void factorize(const MatrixType1& A, MatrixType2& Q, MatrixType3& H);

    /**
     * @brief Computes the Hessenberg factorization of a matrix.
     * @param A The matrix which will be factorized.
     * @param H At the end, an upper Hessenberg matrix with nonnegative
     *        elements on its first subdiagonal.
     * @note The Hessenberg factorization of @c A is
     *       <tt>A = QHQ<sup>t</sup></tt> for an orthogonal matrix @c Q, but
     *       this function only computes @c H.
     */
    template<typename MatrixType1, typename MatrixType2>
    static void factorize(const MatrixType1& A, MatrixType2& H);

}; /* HessenbergFactorization */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename MatrixType1, typename MatrixType2, typename MatrixType3>
void HessenbergFactorization::factorize(const MatrixType1& A,
                                        MatrixType2& Q,
                                        MatrixType3& H)
{
    using T = typename MatrixType1::value_type;
    using size_type = typename MatrixType1::size_type;

    FREEAML_ASSERT(A.is_square() == true);

    const size_type n = A.num_rows();

    /*
     * assign the correct dimensions to Q (at the end, Q will be the transpose
     * of a matrix which is the product of Householder transformations)
     */
    Q = MatrixType2(n, n);

    /* initialize Q as the n × n identity matrix */
    for (size_type i = 0; i < n; ++i)
    {
        Q(i, i) = T{1};
    }

    /* H <-- A */
    H = MatrixType3(n, n, A.flatten());

    /*
     * for each column p of H which must have the portion below its subdiagonal
     * zeroed out
     */
    for (size_type p = 0; p + 1 < n; ++p)
    {
        size_type N = n - p - 1;

        /* let u := H(p+1:n-1, p) */
        Vector<T> u(N);

        for (size_type i = 0; i < N; ++i)
        {
            u[i] = H(i + p + 1, p);
        }

        /*
         * define the Householder transformation for u which projects it onto
         * e_0 = (1, 0, ..., 0)
         */
        Vector<T> v;

        T b, c = 0;

        householder(u, v, b, c, 0);

        /* If b = 0, no Householder transformation is necessary */
        if (b != T{0})
        {
            Matrix<T> B(N, N + 1);
            Matrix<T> C(N, n);

            for (size_type i = 0; i < N; ++i)
            {
                /* B <-- (I - bvv^t)H(p+1:n-1, p:n-1) */
                for (size_type j = 0; j <= N; ++j)
                {
                    B(i, j) = H(i + p + 1, j + p);

                    for (size_type k = 0; k < N; ++k)
                    {
                        B(i, j) -= b * v[i] * v[k] * H(k + p + 1, j + p);
                    }
                }

                /* C <-- (I - bvv^t)Q(0:n-1, p+1:n-1) */
                for (size_type j = 0; j < n; ++j)
                {
                    C(i, j) = Q(j, i + p + 1);

                    for (size_type k = 0; k < N; ++k)
                    {
                        C(i, j) -= b * v[i] * v[k] * Q(j, k + p + 1);
                    }
                }
            }

            for (size_type i = 0; i < N; ++i)
            {
                /* H(p+1:n-1, p:n-1) <-- -c*B */
                for (size_type j = 0; j <= N; ++j)
                {
                    H(i + p + 1, j + p) = -c * B(i, j);
                }

                /* Q(0:n-1, p+1:n-1) <-- -c*C */
                for (size_type j = 0; j < n; ++j)
                {
                    Q(j, i + p + 1) = -c * C(i, j);
                }
            }

            Matrix<T> D(n, N);

            /* D <-- H(0:n-1, p+1:n-1)(I - b*vv^t) */
            for (size_type i = 0; i < n; ++i)
            {
                for (size_type j = 0; j < N; ++j)
                {
                    D(i, j) = H(i, j + p + 1);

                    for (size_type k = 0; k < N; ++k)
                    {
                        D(i, j) -= b * H(i, k + p + 1) * v[k] * v[j];
                    }
                }
            }

            /* H(0:n-1, p+1:n-1) <-- -c*D */
            for (size_type i = 0; i < n; ++i)
            {
                for (size_type j = 0; j < N; ++j)
                {
                    H(i, j + p + 1) = -c * D(i, j);
                }
            }
        }
    }
}

template<typename MatrixType1, typename MatrixType2>
void HessenbergFactorization::factorize(const MatrixType1& A, MatrixType2& H)
{
    using T = typename MatrixType1::value_type;
    using size_type = typename MatrixType1::size_type;

    FREEAML_ASSERT(A.is_square() == true);

    const size_type n = A.num_rows();

    /* H <-- A */
    H = MatrixType2(n, n, A.flatten());

    /*
     * for each column p of H which must have the portion below its subdiagonal
     * zeroed out
     */
    for (size_type p = 0; p + 1 < n; ++p)
    {
        size_type N = n - p - 1;

        /* let u := H(p+1:n-1, p) */
        Vector<T> u(N);

        for (size_type i = 0; i < N; ++i)
        {
            u[i] = H(i + p + 1, p);
        }

        /*
         * define the Householder transformation for u which projects it onto
         * e_0 = (1, 0, ..., 0)
         */
        Vector<T> v;

        T b, c = 0;

        householder(u, v, b, c, 0);

        /* If b = 0, no Householder transformation is necessary */
        if (b != T{0})
        {
            Matrix<T> B(N, N + 1);

            for (size_type i = 0; i < N; ++i)
            {
                /* B <-- (I - bvv^t)H(p+1:n-1, p:n-1) */
                for (size_type j = 0; j <= N; ++j)
                {
                    B(i, j) = H(i + p + 1, j + p);

                    for (size_type k = 0; k < N; ++k)
                    {
                        B(i, j) -= b * v[i] * v[k] * H(k + p + 1, j + p);
                    }
                }
            }

            for (size_type i = 0; i < N; ++i)
            {
                /* H(p+1:n-1, p:n-1) <-- -c*B */
                for (size_type j = 0; j <= N; ++j)
                {
                    H(i + p + 1, j + p) = -c * B(i, j);
                }
            }

            Matrix<T> D(n, N);

            /* D <-- H(0:n-1, p+1:n-1)(I - b*vv^t) */
            for (size_type i = 0; i < n; ++i)
            {
                for (size_type j = 0; j < N; ++j)
                {
                    D(i, j) = H(i, j + p + 1);

                    for (size_type k = 0; k < N; ++k)
                    {
                        D(i, j) -= b * H(i, k + p + 1) * v[k] * v[j];
                    }
                }
            }

            /* H(0:n-1, p+1:n-1) <-- -c*D */
            for (size_type i = 0; i < n; ++i)
            {
                for (size_type j = 0; j < N; ++j)
                {
                    H(i, j + p + 1) = -c * D(i, j);
                }
            }
        }
    }
}

} /* namespace freeaml */
