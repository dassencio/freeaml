#ifndef _freeAML_normal_eqs_h_
#define _freeAML_normal_eqs_h_


namespace aml
{


template < class T >
class normal_eqs
{

public:

	/**
	 * @brief solves a least squares problem using the normal equations
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @param lss a linear system solver class
	 * @return true if the residual tolerance could be achieved within the
	 *         maximum number of iterations allowed, false otherwise
	 */
	template < class _mat, class _vec, class _sys_solver >
	bool solve (const _mat& A, _vec& x, const _vec& b, _sys_solver& lss);

}; /* end of class normal_eqs */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
template < class _mat, class _vec, class _sys_solver >
bool normal_eqs< T >::solve (const _mat&  A,
                                   _vec&  x,
                             const _vec&  b,
                             _sys_solver& lss)
{
	_mat At = A.transpose();

	/* Solve (At*A)*x = At*b */
	return lss.solve(At*A, x, At*b);
}


} /* end of namespace aml */

#endif /* _freeAML_normal_eqs_h_ */
