#pragma once

#include <IterativeLinearSystemSolverBase.h>
#include <rotation.h>

namespace freeaml
{
template<typename T>
class MinimumResidual : public IterativeLinearSystemSolverBase<T>
{
public:
    using BaseSolver = IterativeLinearSystemSolverBase<T>;
    using size_type = typename BaseSolver::size_type;

    /**
     * @brief Constructs a linear system solver which uses the minimum residual
     *        method (MINRES).
     * @param max_iterations The maximum number of iterations allowed.
     * @param residual_tolerance The residual tolerance.
     */
    MinimumResidual(size_type max_iterations, const T& residual_tolerance);

    /**
     * @brief Solves a linear system using the minimum residual method (MINRES).
     * @param A The linear system matrix.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear system.
     * @return @c true if the residual tolerance could be achieved within the
     *         maximum number of iterations allowed, @c false otherwise.
     */
    template<typename MatrixType, typename VectorType>
    bool solve(const MatrixType& A, VectorType& x, const VectorType& b);

}; /* MinimumResidual<T> */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename T>
MinimumResidual<T>::MinimumResidual(const size_type max_iterations,
                                    const T& residual_tolerance)
    : BaseSolver(max_iterations, residual_tolerance)
{
    /* nothing needs to be done here */
}

template<typename T>
template<typename MatrixType, typename VectorType>
bool MinimumResidual<T>::solve(const MatrixType& A,
                               VectorType& x,
                               const VectorType& b)
{
    FREEAML_ASSERT(BaseSolver::check_dimensions(A, x, b) == true);

    size_type n = A.num_rows();

    /* in this algorithm, the initial guess is the zero vector */
    std::fill(x.begin(), x.end(), T{0});

    /*  q_{k-1}, q_k and q_{k+1}, k is the iteration number */
    VectorType q1(n, T{0});
    VectorType q2(n, T{0});
    VectorType q3(n, T{0});

    /* m_{k-2}, m_{k-1} and m_k */
    VectorType m1(n, T{0});
    VectorType m2(n, T{0});
    VectorType m3(n, T{0});

    /*
     * alpha_k, beta_k, beta_{k+1} (these represent the only nonzero elements
     * of the last column of H_k)
     */
    T a1 = T{0};
    T b1 = T{0};
    T b2 = T{0};

    /*
     * technically, we have here k = -1, q_{k+1} = q_0 = b is the first column
     * of the Kryvlov subspace basis (matrix) Q
     */
    q3 = b;

    b1 = q3.l2_norm();

    (*this).num_iterations_ = 0;

    /* if the right-hand side b is zero, the solution is then x = 0 */
    if (b1 == T{0})
    {
        return true;
    }

    /* normalize q3 */
    q3 /= b1;

    /* t_k, the k-th component of the vector t of the MINRES method */
    T t1 = b1;

    /* Givens rotation parameters which define G(k-2,k-1) */
    T c1 = T{0};
    T s1 = T{0};

    /* Givens rotation parameters which define G(k-1,k) */
    T c2 = T{0};
    T s2 = T{0};

    /* Givens rotation parameters which define G(k,k+1) */
    T c3 = T{0};
    T s3 = T{0};

    while ((*this).num_iterations_ < (*this).max_iterations())
    {
        /* q_{k-1} <-- q_k and q_k <-- q_{k+1} */
        q1 = q2;
        q2 = q3;

        /* q_{k+1} <-- A*q_k */
        q3 = A * q2;

        /* alpha_k <-- (q_k)^t*q_{k+1} */
        a1 = q2 * q3;

        /* q_{k+1} <-- q_{k+1} - alpha_k*q_k - beta_k*q_{k-1} */
        q3 -= a1 * q2 + b1 * q1;

        /*
         * beta_{k+1} <-- norm(q_{k+1})_2; last column of H_k is beta_k,
         * alpha_k, beta_{k+1}
         */
        b2 = q3.l2_norm();

        q3 /= b2;

        /* let (epsilon_k,delta_k,gamma_k) = (0,beta_k,alpha_k) */
        T ek = T{0};
        T dk = b1;
        T gk = a1;

        /* make a copy zeta_k of beta_{k+1} */
        T zk = b2;

        m3 = q2;

        if ((*this).num_iterations_ > 1)
        {
            /*
             * rotate (epsilon_k,delta_k) using G(k-2,k-1) to obtain epsilon_k;
             * delta_k is temporary
             */
            givens_rotation(c1, s1, ek, dk);

            m3 -= ek * m1;
        }
        if ((*this).num_iterations_ > 0)
        {
            /*
             * rotate (delta_k,gamma_k) using G(k-1,k) to obtain delta_k;
             * gamma_k is temporary
             */
            givens_rotation(c2, s2, dk, gk);

            m3 -= dk * m2;
        }

        /*
         * get the Givens rotation matrix parameters to make the term T(k+1,k)
         * equal to zero
         */
        givens(gk, b2, c3, s3);

        /*
         * rotate (gamma_k,zeta_k) using G(k,k+1) to obtain gamma_k and zero out
         * beta_{k+1} on T
         */
        givens_rotation(c3, s3, gk, zk);

        m3 /= gk;

        /*
         * t_{k+1}, the (k+1)-th component of the vector t of the MINRES
         * method
         */
        T t2 = T{0};

        /*
         * apply the Givens rotation G(k,k+1) to the last two nonzero elements
         * of the right-hand side: t_{k+1} <-- t_k*s_k and t_k <-- t_k*c_k
         */
        givens_rotation(c3, s3, t1, t2);

        x += t1 * m3;

        /* t_k <-- t_{k+1} */
        t1 = t2;

        /* m_{k-2} <-- m_{k-1} and m_{k-1} <-- m_k */
        m1 = m2;
        m2 = m3;

        c1 = c2;
        s1 = s2;

        c2 = c3;
        s2 = s3;

        /* beta_k <-- beta_{k+1} */
        b1 = b2;

        ++(*this).num_iterations_;

        /* if the residual is within the specified tolerance, stop */
        if (std::abs(t1) <= (*this).residual_tolerance())
        {
            return true;
        }
    }

    return false;
}

} /* namespace freeaml */
