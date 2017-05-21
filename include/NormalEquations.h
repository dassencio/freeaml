#pragma once

namespace freeaml
{
class NormalEquations
{
public:
    /**
     * @brief Solves a linear least-squares problem using its normal equations.
     * @param A The matrix of the linear least-squares problem.
     * @param x The vector on which the solution will be written.
     * @param b The right-hand side of the linear least-squares problem.
     * @param lss A linear system solver.
     * @return @c true if a solution could be found, @c false otherwise.
     */
    template<typename MatrixType, typename VectorType, typename LinSysSolver>
    static bool solve(const MatrixType& A,
                      VectorType& x,
                      const VectorType& b,
                      LinSysSolver& lss);

}; /* NormalEquations */

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename MatrixType, typename VectorType, typename LinSysSolver>
bool NormalEquations::solve(const MatrixType& A,
                            VectorType& x,
                            const VectorType& b,
                            LinSysSolver& lss)
{
    MatrixType At = A.transpose();

    /* solve (A^t*A)*x = A^t*b */
    return lss.solve(At * A, x, At * b);
}

} /* namespace freeaml */
