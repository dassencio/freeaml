#pragma once

#include <IterativeLinearSystemSolverBase.h>

namespace freeaml
{
template<typename T>
class Jacobi : public IterativeLinearSystemSolverBase<T>
{
public:
    using BaseSolver = IterativeLinearSystemSolverBase<T>;
    using size_type = typename BaseSolver::size_type;

    /**
     * @brief Constructs a linear system solver which uses the Jacobi method.
     * @param max_iterations The maximum number of iterations allowed.
     * @param residual_tolerance The residual tolerance.
     */
    Jacobi(size_type max_iterations, const T& residual_tolerance);

    /**
     * @brief Solves a linear system using the Jacobi method.
     * @param A The linear system matrix.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear system.
     * @return @c true if the residual tolerance could be achieved within the
     *         maximum number of iterations allowed, @c false otherwise.
     */
    template<typename MatrixType, typename VectorType>
    bool solve(const MatrixType& A, VectorType& x, const VectorType& b);

    /**
     * @brief Performs a fixed number of Jacobi iterations for a linear system.
     * @param A The linear system matrix.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear system.
     * @param n The number of iterations that will be performed.
     * @note This function simply performs @c n Jacobi iterations without taking
     *       the residual tolerance into account.
     */
    template<typename MatrixType, typename VectorType>
    static void iterate(const MatrixType& A,
                        VectorType& x,
                        const VectorType& b,
                        size_type n);

}; /* class Jacobi<T> */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename T>
Jacobi<T>::Jacobi(const size_type max_iterations, const T& residual_tolerance)
    : BaseSolver(max_iterations, residual_tolerance)
{
    /* nothing needs to be done here */
}

template<typename T>
template<typename MatrixType, typename VectorType>
bool Jacobi<T>::solve(const MatrixType& A, VectorType& x, const VectorType& b)
{
    FREEAML_ASSERT(BaseSolver::check_dimensions(A, x, b));

    /* r = residual vector */
    VectorType r = b - A * x;

    (*this).num_iterations_ = 0;

    while ((*this).num_iterations_ < (*this).max_iterations())
    {
        for (size_type i = 0; i < x.size(); ++i)
        {
            FREEAML_ASSERT(A(i, i) != T{0});

            x[i] += r[i] / A(i, i);
        }

        r = b - A * x;

        ++(*this).num_iterations_;

        /* if the residual is within the maximum tolerance, stop */
        if (r.l2_norm() <= (*this).residual_tolerance())
        {
            return true;
        }
    }

    return false;
}

template<typename T>
template<typename MatrixType, typename VectorType>
void Jacobi<T>::iterate(const MatrixType& A,
                        VectorType& x,
                        const VectorType& b,
                        const size_type n)
{
    FREEAML_ASSERT(BaseSolver::check_dimensions(A, x, b));

    /* r = residual vector */
    VectorType r = b - A * x;

    size_type num_iterations = 0;

    while (num_iterations < n)
    {
        for (size_type i = 0; i < x.size(); ++i)
        {
            x[i] += r[i] / A(i, i);
        }

        ++num_iterations;

        r = b - A * x;
    }
}

} /* namespace freeaml */
