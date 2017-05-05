#pragma once

#include <IterativeLinearSystemSolverBase.h>

namespace freeaml
{
template<typename T>
class BiconjugateGradientStabilized : public IterativeLinearSystemSolverBase<T>
{
public:
    using BaseSolver = IterativeLinearSystemSolverBase<T>;
    using size_type = typename BaseSolver::size_type;

    /**
     * @brief Constructs a linear system solver which uses the biconjugate
     *        gradient stabilized method (BiCGSTAB).
     * @param max_iterations The maximum number of iterations allowed.
     * @param residual_tolerance The maximum residual tolerance allowed.
     */
    BiconjugateGradientStabilized(size_type max_iterations,
                                  const T& residual_tolerance);

    /**
     * @brief Solves a linear system using the biconjugate gradient stabilized
     *        method (BiCGSTAB).
     * @param A The linear system matrix.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear system.
     * @return @c true if the residual tolerance could be achieved within the
     *         maximum number of iterations allowed, @c false otherwise.
     */
    template<typename MatrixType, typename VectorType>
    bool solve(const MatrixType& A, VectorType& x, const VectorType& b);

}; /* class BiconjugateGradientStabilized<T> */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename T>
BiconjugateGradientStabilized<T>::BiconjugateGradientStabilized(
    const size_type max_iterations, const T& residual_tolerance)
    : BaseSolver(max_iterations, residual_tolerance)
{
    /* nothing needs to  be done here */
}

template<typename T>
template<typename MatrixType, typename VectorType>
bool BiconjugateGradientStabilized<T>::solve(const MatrixType& A,
                                             VectorType& x,
                                             const VectorType& b)
{
    FREEAML_ASSERT(BaseSolver::check_dimensions(A, x, b));

    /* r = residual vector */
    VectorType r = b - A * x;

    VectorType u = r;

    T rho{1};
    T alpha{1};
    T omega{1};

    VectorType p(b.size(), T{0});
    VectorType v(b.size(), T{0});

    (*this).num_iterations_ = 0;

    while ((*this).num_iterations_ < (*this).max_iterations())
    {
        T old_rho = rho;

        rho = u * r;

        T beta = (rho * alpha) / (old_rho * omega);

        p = r + beta * (p - omega * v);
        v = A * p;

        alpha = rho / (u * v);

        VectorType s = r - alpha * v;
        VectorType t = A * s;

        omega = (t * s) / (t * t);

        x += alpha * p + omega * s;

        r = s - omega * t;

        ++(*this).num_iterations_;

        /* if the residual is within the maximum tolerance, stop */
        if (r.l2_norm() <= (*this).residual_tolerance())
        {
            return true;
        }
    }

    return false;
}

} /* namespace freeaml */
