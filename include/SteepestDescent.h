#pragma once

#include <IterativeLinearSystemSolverBase.h>
#include <cmath>

namespace freeaml
{
template<typename T>
class SteepestDescent : public IterativeLinearSystemSolverBase<T>
{
public:
    using BaseSolver = IterativeLinearSystemSolverBase<T>;
    using size_type = typename BaseSolver::size_type;

    /**
     * @brief Constructs a linear system solver which uses the steepest descent
     *        method.
     * @param max_iterations The maximum number of iterations allowed.
     * @param residual_tolerance The residual tolerance.
     */
    SteepestDescent(size_type max_iterations, const T& residual_tolerance);

    /**
     * @brief Solves a linear system using the steepest descent method.
     * @param A The linear system matrix.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear system.
     * @return @c true if the residual tolerance could be achieved within the
     *         maximum number of iterations allowed, @c false otherwise.
     */
    template<typename MatrixType, typename VectorType>
    bool solve(const MatrixType& A, VectorType& x, const VectorType& b);

}; /* class SteepestDescent<T> */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename T>
SteepestDescent<T>::SteepestDescent(const size_type max_iterations,
                                    const T& residual_tolerance)
    : BaseSolver(max_iterations, residual_tolerance)
{
    /* nothing needs to be done here */
}

template<typename T>
template<typename MatrixType, typename VectorType>
bool SteepestDescent<T>::solve(const MatrixType& A,
                               VectorType& x,
                               const VectorType& b)
{
    FREEAML_ASSERT(BaseSolver::check_dimensions(A, x, b) == true);

    /* r = residual vector */
    VectorType r = b - A * x;

    T rr = r * r;

    (*this).num_iterations_ = 0;

    while ((*this).num_iterations_ < (*this).max_iterations())
    {
        VectorType Ar = A * r;

        T alpha = rr / (r * Ar);

        x += alpha * r;
        r -= alpha * Ar;

        rr = r * r;

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
