#pragma once

#include <debug.h>

namespace freeaml
{
/**
 * @brief Base class defining the properties required by all iterative linear
 *        system solvers.
 */
template<typename T>
class IterativeLinearSystemSolverBase
{
public:
    using size_type = size_t;

    /**
     * @brief Constructs an iterative linear system solver.
     * @param max_iterations The maximum number of iterations allowed.
     * @param residual_tolerance The residual tolerance.
     */
    IterativeLinearSystemSolverBase(size_type max_iterations,
                                    const T& residual_tolerance);

    /**
     * @brief Gets the number of iterations performed in the last solve.
     * @return The number of iterations performed in the last solve.
     */
    size_type num_iterations() const;

    /**
     * @brief Gets the maximum number of iterations allowed in one solve.
     * @return The maximum number of iterations allowed in one solve.
     */
    size_type max_iterations() const;

    /**
     * @brief Gets the residual tolerance limit for the solver.
     * @return The residual tolerance limit for the solver.
     */
    const T& residual_tolerance() const;

    /**
     * @brief Checks if the dimensions of the matrix and vectors which define a
     *        linear system are correct.
     * @param A The linear system matrix.
     * @param x The solution vector.
     * @param b The right-hand side of the linear system.
     * @return @c true if @c A is square and has the same number of rows as @c x
     *         and @c b, @c false otherwise.
     */
    template<typename MatrixType, typename VectorType>
    static bool check_dimensions(const MatrixType& A,
                                 const VectorType& x,
                                 const VectorType& b);

protected:
    size_type num_iterations_;       /* number of iterations performed */
    const size_type max_iterations_; /* maximum number of iterations allowed */
    const T residual_tolerance_;     /* residual tolerance for the solver */

}; /* class IterativeLinearSystemSolverBase<T> */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename T>
IterativeLinearSystemSolverBase<T>::IterativeLinearSystemSolverBase(
    const size_type max_iterations, const T& residual_tolerance)
    : num_iterations_{0},
      max_iterations_(max_iterations),
      residual_tolerance_(residual_tolerance)
{
    FREEAML_ASSERT(max_iterations > 0);
    FREEAML_ASSERT(residual_tolerance > T{0});
}

template<typename T>
typename IterativeLinearSystemSolverBase<T>::size_type
IterativeLinearSystemSolverBase<T>::num_iterations() const
{
    return num_iterations_;
}

template<typename T>
typename IterativeLinearSystemSolverBase<T>::size_type
IterativeLinearSystemSolverBase<T>::max_iterations() const
{
    return max_iterations_;
}

template<typename T>
const T& IterativeLinearSystemSolverBase<T>::residual_tolerance() const
{
    return residual_tolerance_;
}

template<typename T>
template<typename MatrixType, typename VectorType>
bool IterativeLinearSystemSolverBase<T>::check_dimensions(const MatrixType& A,
                                                          const VectorType& x,
                                                          const VectorType& b)
{
    return A.is_square() && A.num_rows() == x.size() && x.size() == b.size();
}

} /* namespace freeaml */
