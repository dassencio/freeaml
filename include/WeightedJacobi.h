#pragma once

#include <IterativeLinearSystemSolverBase.h>

namespace freeaml
{
template<typename T>
class WeightedJacobi : public IterativeLinearSystemSolverBase<T>
{
public:
    using BaseSolver = IterativeLinearSystemSolverBase<T>;
    using size_type = typename BaseSolver::size_type;

    /**
     * @brief Constructs a linear system solver which uses the weighted-Jacobi
     *        method.
     * @param max_iterations The maximum number of iterations allowed.
     * @param residual_tolerance The residual tolerance.
     * @param omega The iteration parameter (weight).
     */
    WeightedJacobi(size_type max_iterations,
                   const T& residual_tolerance,
                   const T& omega);

    /**
     * @brief Solves a linear system using the weighted-Jacobi method.
     * @param A The linear system matrix.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear system.
     * @return @c true if the residual tolerance could be achieved within the
     *         maximum number of iterations allowed, @c false otherwise.
     */
    template<typename MatrixType, typename VectorType>
    bool solve(const MatrixType& A, VectorType& x, const VectorType& b);

    /**
     * @brief Performs a fixed number of weighted-Jacobi iterations for a
     *        linear system.
     * @param A The linear system matrix.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear system.
     * @param n The number of iterations that will be performed.
     * @param omega The iteration parameter (weight).
     * @note This function simply performs @c n weighted-Jacobi iterations
     *       without taking the residual tolerance into account.
     */
    template<typename MatrixType, typename VectorType>
    static void iterate(const MatrixType& A,
                        VectorType& x,
                        const VectorType& b,
                        size_type n,
                        const T& omega);

    /**
     * @brief Gets the iteration parameter (weight) set for the solver.
     * @return The iteration parameter (weight) set for the solver.
     */
    const T& omega() const;

private:
    const T omega_;

}; /* class WeightedJacobi<T> */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename T>
WeightedJacobi<T>::WeightedJacobi(const size_type max_iterations,
                                  const T& residual_tolerance,
                                  const T& omega)
    : BaseSolver(max_iterations, residual_tolerance), omega_(omega)
{
    FREEAML_ASSERT(omega > T{0});
}

template<typename T>
template<typename MatrixType, typename VectorType>
bool WeightedJacobi<T>::solve(const MatrixType& A,
                              VectorType& x,
                              const VectorType& b)
{
    FREEAML_ASSERT(BaseSolver::check_dimensions(A, x, b) == true);

    /* r = residual vector */
    VectorType r = b - A * x;

    (*this).num_iterations_ = 0;

    while ((*this).num_iterations_ < (*this).max_iterations())
    {
        for (size_type i = 0; i < x.size(); ++i)
        {
            FREEAML_ASSERT(A(i, i) != T{0});

            x[i] += omega() * r[i] / A(i, i);
        }

        r = b - A * x;

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
template<typename MatrixType, typename VectorType>
void WeightedJacobi<T>::iterate(const MatrixType& A,
                                VectorType& x,
                                const VectorType& b,
                                const size_type n,
                                const T& omega)
{
    FREEAML_ASSERT(BaseSolver::check_dimensions(A, x, b) == true);

    /* r = residual vector */
    VectorType r = b - A * x;

    for (size_type num_iterations = 0; num_iterations < n; ++num_iterations)
    {
        for (size_type i = 0; i < x.size(); ++i)
        {
            FREEAML_ASSERT(A(i, i) != T{0});

            x[i] += omega * r[i] / A(i, i);
        }

        r = b - A * x;
    }
}

template<typename T>
const T& WeightedJacobi<T>::omega() const
{
    return omega_;
}

} /* namespace freeaml */
