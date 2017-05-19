#pragma once

#include <IterativeLinearSystemSolverBase.h>
#include <Matrix.h>
#include <rotation.h>

namespace freeaml
{
template<typename T>
class GeneralizedMinimumResidual : public IterativeLinearSystemSolverBase<T>
{
public:
    using BaseSolver = IterativeLinearSystemSolverBase<T>;
    using size_type = typename BaseSolver::size_type;

    /**
     * @brief Constructs a linear system solver which uses the generalized
     *        minimum residual method (GMRES).
     * @param max_iterations The maximum number of iterations allowed.
     * @param residual_tolerance The residual tolerance.
     * @note This GMRES implementation is based on Givens rotations.
     */
    GeneralizedMinimumResidual(size_type max_iterations,
                               const T& residual_tolerance);

    /**
     * @brief Solves a linear system using the generalized minimum residual
     *        method (GMRES).
     * @param A The linear system matrix.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear system.
     * @return @c true if the residual tolerance could be achieved within the
     *         maximum number of iterations allowed, @c false otherwise.
     */
    template<typename MatrixType, typename VectorType>
    bool solve(const MatrixType& A, VectorType& x, const VectorType& b);

}; /* class GeneralizedMinimumResidual<T> */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename T>
GeneralizedMinimumResidual<T>::GeneralizedMinimumResidual(
    const size_type max_iterations, const T& residual_tolerance)
    : BaseSolver(max_iterations, residual_tolerance)
{
    /* nothing needs to be done here */
}

template<typename T>
template<typename MatrixType, typename VectorType>
bool GeneralizedMinimumResidual<T>::solve(const MatrixType& A,
                                          VectorType& x,
                                          const VectorType& b)
{
    FREEAML_ASSERT(BaseSolver::check_dimensions(A, x, b) == true);

    const size_type n = A.num_rows();

    Matrix<T> Q(1, n);
    Matrix<T> H(1, n);

    const T norm_b = b.l2_norm();

    /* build q_0, the basis of the Krylov subspace K^1 = {b} */
    for (size_type j = 0; j < n; ++j)
    {
        Q(0, j) = b[j] / norm_b;
    }

    VectorType f(n);

    /* f = (norm_b, 0, ... 0) */
    f[0] = norm_b;

    /*
     * c and s will hold the parameters needed to define the successive Givens
     * rotations applied (not all elements of these vectors will be used if
     * convergence happens in less than n iterations)
     */
    VectorType c(n);
    VectorType s(n);

    (*this).num_iterations_ = 0;

    while ((*this).num_iterations_ < (*this).max_iterations())
    {
        size_type k = (*this).num_iterations_;

        Q.resize(k + 2, n);
        H.resize(k + 2, n);

        /* build q_k as the k-th row of Q */
        VectorType q_k(n);

        for (size_type j = 0; j < n; ++j)
        {
            q_k[j] = Q(k, j);
        }

        /*
         * z will be used to compute q_{k+1} if we have not yet built the
         * largest Krylov subspace (the one with n dimensions)
         */
        VectorType z = A * q_k;

        /*
         * project out the components of z parallel to {q_0,q_1,...,q_k} and
         * build h_k (the k-th column of H)
         */
        for (size_type i = 0; i <= k; ++i)
        {
            for (size_type j = 0; j < n; ++j)
            {
                H(i, k) += Q(i, j) * z[j];
            }

            for (size_type j = 0; j < n; ++j)
            {
                z[j] -= H(i, k) * Q(i, j);
            }
        }

        /*
         * apply the previous Givens rotations that kept H upper triangular (H
         * here being the theoretical H)
         */
        for (size_type i = 0; i < k; ++i)
        {
            givens_rotation(c[i], s[i], H(i, k), H(i + 1, k));
        }

        /* if one more Givens rotation must be applied */
        if (k + 1 < n)
        {
            T norm_z = z.l2_norm();

            for (size_type j = 0; j < n; ++j)
            {
                Q(k + 1, j) = z[j] / norm_z;
            }

            H(k + 1, k) = norm_z;

            /*
             * compute the parameters for the Givens rotation matrix G(k,k+1)
             * that zeros out H(k+1,k)
             */
            givens(H(k, k), H(k + 1, k), c[k], s[k]);

            /* keep H (the theoretical H) upper triangular */
            givens_rotation(c[k], s[k], H(k, k), H(k + 1, k));

            /* apply also the Givens rotation to the right-hand side f */
            givens_rotation(c[k], s[k], f[k], f[k + 1]);
        }

        /* if H(k,k) = 0, A is not invertible */
        if (H(k, k) == T{0})
        {
            return false;
        }

        ++(*this).num_iterations_;

        /*
         * if we can no longer iterate or if the residual is within the
         * specified tolerance
         */
        if (k + 1 == n || std::abs(f[k + 1]) <= (*this).residual_tolerance())
        {
            VectorType lambda(k + 1);

            /* solve H_{0:k,0:k}*lambda = f */
            size_type i = k + 1;

            while (i-- > 0)
            {
                T y = T{0};

                for (size_type j = i + 1; j <= k; ++j)
                {
                    y += H(i, j) * lambda[j];
                }

                lambda[i] = (f[i] - y) / H(i, i);
            }

            /*
             * build the problem solution x = Q*lambda (in this implementation,
             * Q is the theoretical Q^t!)
             */
            std::fill(x.begin(), x.end(), T{0});

            for (size_type j = 0; j < n; ++j)
            {
                for (size_type i = 0; i <= k; ++i)
                {
                    x[j] += lambda[i] * Q(i, j);
                }
            }

            return (b - A * x).l2_norm() <= (*this).residual_tolerance();
        }
    }

    return false;
}

} /* namespace aml */
