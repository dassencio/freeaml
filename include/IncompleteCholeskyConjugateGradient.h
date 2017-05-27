#pragma once

#include <IterativeLinearSystemSolverBase.h>
#include <SparseMatrix.h>
#include <cmath>

namespace freeaml
{
template<typename T>
class IncompleteCholeskyConjugateGradient
    : public IterativeLinearSystemSolverBase<T>
{
public:
    using BaseSolver = IterativeLinearSystemSolverBase<T>;
    using size_type = typename BaseSolver::size_type;

    /**
     * @brief Constructs a linear system solver which uses the conjugate
     *        gradient method (CG) with incomplete-Cholesky preconditioning.
     * @param max_iterations The maximum number of iterations allowed.
     * @param residual_tolerance The residual tolerance.
     */
    IncompleteCholeskyConjugateGradient(size_type max_iterations,
                                        const T& residual_tolerance);

    /**
     * @brief Solves a linear system using the conjugate gradient method (CG)
     *        with incomplete-Cholesky preconditioning.
     * @param A The linear system matrix.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear system.
     * @return @c true if the residual tolerance could be achieved within the
     *         maximum number of iterations allowed, @c false otherwise.
     * @note This method should only be used if @c A is sparse.
     */
    template<typename MatrixType, typename VectorType>
    bool solve(const MatrixType& A, VectorType& x, const VectorType& b);

private:
    /**
     * @brief Builds the (sparse) lower triangular matrix @c K which defines the
     *        incomplete Cholesky factorization of a matrix @c A, i.e., the
     *        matrix @c K such that <tt>KK<sup>t</sup></tt> approximates @c A.
     * @param A A matrix.
     * @param K The matrix on which the incomplete Cholesky factorization matrix
     *        @c K will be written (@c A is assumed to be sparse).
     * @return true if @c K could be built successfully, @c false otherwise.
     */
    template<typename MatrixType>
    bool build_preconditioner(const MatrixType& A, SparseMatrix<T>& K) const;

    /**
     * @brief Solves the incomplete-Cholesky preconditioner equation
     *        <tt>KK<sup>t</sup>z = r</tt>.
     * @param K A sparse lower triangular matrix.
     * @param z The vector on which the solution will be written.
     * @param r The right-hand side vector (in our case, a residual vector).
     *
     */
    template<typename VectorType>
    void solve_preconditioner_equation(const SparseMatrix<T>& K,
                                       VectorType& z,
                                       const VectorType& r) const;

}; /* class IncompleteCholeskyConjugateGradient<T> */

/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template<typename T>
IncompleteCholeskyConjugateGradient<T>::IncompleteCholeskyConjugateGradient(
    const size_type max_iterations, const T& residual_tolerance)
    : BaseSolver(max_iterations, residual_tolerance)
{
    /* nothing needs to be done here */
}

template<typename T>
template<typename MatrixType, typename VectorType>
bool IncompleteCholeskyConjugateGradient<T>::solve(const MatrixType& A,
                                                   VectorType& x,
                                                   const VectorType& b)
{
    FREEAML_ASSERT(BaseSolver::check_dimensions(A, x, b) == true);

    (*this).num_iterations_ = 0;

    const size_type n = A.num_rows();

    SparseMatrix<T> K(n, n);

    /* the incomplete Cholesky preconditioner is KK^t */
    if (build_preconditioner(A, K) == false)
    {
        return false;
    }

    /* r = residual vector */
    VectorType r = b - A * x;

    VectorType z(n);

    /* z = (KK^t)^(-1)r */
    solve_preconditioner_equation(K, z, r);

    VectorType p = z;
    VectorType q = A * z;

    while ((*this).num_iterations_ < (*this).max_iterations())
    {
        T gamma = r * z;
        T alpha = gamma / (p * q);

        x += alpha * p;
        r -= alpha * q;

        /* z = (KK^t)^(-1)r */
        solve_preconditioner_equation(K, z, r);

        T beta = (r * z) / gamma;

        p = z + (beta * p);
        q = A * z + (beta * q);

        ++(*this).num_iterations_;

        /* if the residual is within the specified tolerance, stop */
        if (r.l2_norm() <= (*this).residual_tolerance())
        {
            return true;
        }
    }

    return false;
}

template<typename T>
template<typename MatrixType>
bool IncompleteCholeskyConjugateGradient<T>::build_preconditioner(
    const MatrixType& A, SparseMatrix<T>& K) const
{
    const size_type n = A.num_rows();

    for (size_type j = 0; j < n; ++j)
    {
        if (A(j, j) != T{0})
        {
            K(j, j) = A(j, j);

            for (const auto& element : K.row(j))
            {
                if (element.first >= j)
                {
                    break;
                }

                K(j, j) -= element.second * element.second;
            }
        }

        FREEAML_ASSERT(K(j, j) >= T{0});

        K(j, j) = std::sqrt(K(j, j));

        /* if K(j,j) = 0, then A is not invertible */
        if (K(j, j) == T{0})
        {
            return false;
        }

        /* const reference to K (for safe element access) */
        const SparseMatrix<T>& cK = K;

        for (size_type i = j + 1; i < n; ++i)
        {
            if (A(i, j) != T{0})
            {
                K(i, j) = A(i, j);

                for (const auto& element : K.row(i))
                {
                    if (element.first >= j)
                    {
                        break;
                    }

                    K(i, j) -= element.second * cK(j, element.first);
                }

                K(i, j) /= K(j, j);
            }
        }
    }

    return true;
}

template<typename T>
template<typename VectorType>
void IncompleteCholeskyConjugateGradient<T>::solve_preconditioner_equation(
    const SparseMatrix<T>& K, VectorType& z, const VectorType& r) const
{
    /* K is an n × n matrix */
    const size_type n = K.num_rows();

    std::fill(z.begin(), z.end(), T{0});

    VectorType w(n);

    const SparseMatrix<T> Kt = K.transpose();

    /* perform forward substitution to solve Kw = b, with w = (K^t)z */
    for (size_type i = 0; i < n; ++i)
    {
        T y = T{0};

        for (const auto& element : K.row(i))
        {
            if (element.first >= i)
            {
                break;
            }

            y += element.second * w[element.first];
        }

        w[i] = (r[i] - y) / K(i, i);
    }

    /* perform backward substitution to solve (K^t)z = w */
    size_type i = n;

    while (i-- > 0)
    {
        T y = T{0};

        for (const auto& element : Kt.row(i))
        {
            y += element.second * z[element.first];
        }

        z[i] = (w[i] - y) / Kt(i, i);
    }
}

} /* namespace freeaml */
