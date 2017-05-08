#pragma once

#include <IterativeLinearSystemSolverBase.h>

namespace freeaml
{
template<typename T>
class SuccessiveOverRelaxation : public IterativeLinearSystemSolverBase<T>
{
public:
    using BaseSolver = IterativeLinearSystemSolverBase<T>;
    using size_type = typename BaseSolver::size_type;

    /**
     * @brief Constructs a linear system solver which uses the successive
     *        over-relaxation method (SOR).
     * @param max_iterations The maximum number of iterations allowed.
     * @param residual_tolerance The residual tolerance.
     * @param omega The relaxation factor.
     */
    SuccessiveOverRelaxation(size_type max_iterations,
                             const T& residual_tolerance,
                             const T& omega);

    /**
     * @brief Solves a linear system using the successive over-relaxation method
     *        (SOR).
     * @param A The linear system matrix.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear system.
     * @return @c true if the residual tolerance could be achieved within the
     *         maximum number of iterations allowed, @c false otherwise.
     */
    template<typename MatrixType, typename VectorType>
    bool solve(const MatrixType& A, VectorType& x, const VectorType& b);

    /**
     * @brief Gets the relaxation factor set for the solver.
     * @return The relaxation factor set for the solver.
     */
    const T& omega() const;

private:
    const T omega_;

}; /* SuccessiveOverRelaxation<T> */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename T>
SuccessiveOverRelaxation<T>::SuccessiveOverRelaxation(
    const size_type max_iterations, const T& residual_tolerance, const T& omega)
    : BaseSolver(max_iterations, residual_tolerance), omega_(omega)
{
    FREEAML_ASSERT(omega > T{0});
}

template<typename T>
template<typename MatrixType, typename VectorType>
bool SuccessiveOverRelaxation<T>::solve(const MatrixType& A,
                                        VectorType& x,
                                        const VectorType& b)
{
    FREEAML_ASSERT(BaseSolver::check_dimensions(A, x, b));

    VectorType x0 = x;

    (*this).num_iterations_ = 0;

    while ((*this).num_iterations_ < (*this).max_iterations())
    {
        for (size_type i = 0; i < x.size(); ++i)
        {
            /* if A(i,i) = 0, the SOR method cannot be used */
            if (A(i, i) == T{0})
            {
                return false;
            }

            T c = b[i];

            for (size_type j = 0; j < i; ++j)
            {
                c -= A(i, j) * x[j];
            }

            for (size_type j = i + 1; j < A.num_rows(); ++j)
            {
                c -= A(i, j) * x0[j];
            }

            x[i] = (T{1} - omega()) * x0[i] + omega() * (c / A(i, i));
        }

        ++(*this).num_iterations_;

        /* if the residual is within the maximum tolerance, stop */
        if ((b - A * x).l2_norm() <= (*this).residual_tolerance())
        {
            return true;
        }

        x0 = x;
    }

    return false;
}

template<typename T>
const T& SuccessiveOverRelaxation<T>::omega() const
{
    return omega_;
}

} /* namespace aml */
