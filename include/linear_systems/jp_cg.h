#ifndef _freeAML_jp_cg_h_
#define _freeAML_jp_cg_h_


#include "cg.h"

#include <data_structures/sparse_mat.h>


namespace aml
{


template < class T >
class jp_cg: public cg< T >
{

public:

	/***********************************************************************
	 *
	 *	CONSTRUCTORS
	 *
	 **********************************************************************/

	/**
	 * @brief default constructor
	 * @param max_iter maximum number of iterations allowed in one solve
	 * @param residual_tol maximum residual tolerance allowed
	 */
	jp_cg (const size_t max_iter, const T& residual_tol);


	/***********************************************************************
	 *
	 *	SOLVE FUNCTION
	 *
	 **********************************************************************/

	/**
	 * @brief solves a linear system using the Jacobi preconditioned
	 *        conjugate gradient method
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @return true if the residual tolerance could be achieved within the
	 *         maximum number of iterations allowed, false otherwise
	 * @note a Jacobi preconditioner can only be used if all diagonal
	 *       elements of A are nonzero
	 */
	template < class _mat, class _vec >
	bool solve (const _mat& A, _vec& x, const _vec& b);

}; /* end of class jp_cg */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
jp_cg< T >::jp_cg(const size_t max_iter, const T& residual_tol):

        cg< T >(max_iter, residual_tol)
{
	/* nothing needs to  be done here */
}


template < class T >
template < class _mat, class _vec >
bool jp_cg< T >::solve (const _mat& A, _vec& x, const _vec& b)
{
	/* D is the Jacobi preconditioner matrix */
	sparse_mat< T > D(A.num_rows(), A.num_cols());

	cg< T >::set_num_iterations(0);

	/* build D */
	for (size_t i = 0; i < D.num_rows(); i++)
	{
		/* if A(i,i) = 0, the Jacobi preconditioner cannot be built */
		if (A(i,i) == (T) 0)
		{
			return false;
		}
		else
		{
			D(i,i) = (T) 1 / std::sqrt(std::abs(A(i,i)));
		}
	}

	/* solve (DAD)*x = D*b (preconditioned linear system); here x is
	 * not the final solution; D*x is! */
	if (!cg< T >::solve(D*A*D, x, D*b))
	{
		return false;
	}

	x = D*x;

	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_jp_cg_h_ */
