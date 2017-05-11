#pragma once

#include <IterativeLinearSystemSolverBase.h>
#include <cmath>

namespace freeaml
{
template<typename T>
class ConjugateGradient : public IterativeLinearSystemSolverBase<T>
{
public:
    using BaseSolver = IterativeLinearSystemSolverBase<T>;
    using size_type = typename BaseSolver::size_type;

    /**
     * @brief Constructs a linear system solver which uses the conjugate
     *        gradient method (CG).
     * @param max_iterations The maximum number of iterations allowed.
     * @param residual_tolerance The residual tolerance.
     */
    ConjugateGradient(size_type max_iterations, const T& residual_tolerance);

    /**
     * @brief Solves a linear system using the conjugate gradient method (CG).
     * @param A The linear system matrix.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear system.
     * @return @c true if the residual tolerance could be achieved within the
     *         maximum number of iterations allowed, @c false otherwise.
     */
    template<typename MatrixType, typename VectorType>
    bool solve(const MatrixType& A, VectorType& x, const VectorType& b);

}; /* class ConjugateGradient<T> */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename T>
ConjugateGradient<T>::ConjugateGradient(const size_type max_iterations,
                                        const T& residual_tolerance)
    : BaseSolver(max_iterations, residual_tolerance)
{
    /* nothing needs to be done here */
}

template<typename T>
template<typename MatrixType, typename VectorType>
bool ConjugateGradient<T>::solve(const MatrixType& A,
                                 VectorType& x,
                                 const VectorType& b)
{
    FREEAML_ASSERT(BaseSolver::check_dimensions(A, x, b) == true);

    /* r = residual vector */
    VectorType r = b - A * x;

    T rr = r * r;

    VectorType p = r;
    VectorType q = A * r;

    (*this).num_iterations_ = 0;

    while ((*this).num_iterations_ < (*this).max_iterations())
    {
        T alpha = rr / (p * q);

        x += alpha * p;
        r -= alpha * q;

        T old_rr = rr;

        rr = r * r;

        T beta = rr / old_rr;

        p = r + (beta * p);
        q = A * r + (beta * q);

        ++(*this).num_iterations_;

        /* if the residual is within the specified tolerance, stop */
        if (std::sqrt(rr) <= (*this).residual_tolerance())
        {
            return true;
        }
    }

    return false;
}

} /* namespace freeaml */
