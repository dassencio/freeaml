#ifndef _freeAML_cg_h_
#define _freeAML_cg_h_


#include <cmath>

#include "generic_iter.h"


namespace aml
{


template < class T >
class cg: public generic_iter< T >
{

public:

	typedef generic_iter< T >	gen_it_type;


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
	cg (const size_t max_iter, const T& residual_tol);


	/***********************************************************************
	 *
	 *	SOLVE FUNCTION
	 *
	 **********************************************************************/

	/**
	 * @brief solves a linear system using the conjugate gradient method
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @return true if the residual tolerance could be achieved within the
	 *         maximum number of iterations allowed, false otherwise
	 */
	template < class _mat, class _vec >
	bool solve (const _mat& A, _vec& x, const _vec& b);

}; /* end of class cg */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
cg< T >::cg(const size_t max_iter, const T& residual_tol):

        gen_it_type(max_iter, residual_tol)
{
	/* nothing needs to  be done here */
}


template < class T >
template < class _mat, class _vec >
bool cg< T >::solve (const _mat& A, _vec& x, const _vec& b)
{
	FREEAML_ASSERT(gen_it_type::check_dimensions(A,x,b));

	/* r = residual vector */
	_vec r = b - A*x;

	T rr = r*r;

	_vec p = r;
	_vec q = A*r;

	size_t num_iter = 0;

	while (num_iter <= gen_it_type::max_num_iterations())
	{
		T alpha = rr / (p*q);

		x += alpha * p;
		r -= alpha * q;

		T old_rr = rr;

		rr = r*r;

		T beta = rr / old_rr;

		p = r + (beta * p);
		q = A*r + (beta * q);

		num_iter++;

		/* if the residual is within the maximum tolerance, stop */
		if (std::sqrt(rr) <= gen_it_type::residual_tolerance())
		{
			gen_it_type::set_num_iterations(num_iter);
			return true;
		}
	}

	gen_it_type::set_num_iterations(num_iter);
	return false;
}


} /* end of namespace aml */

#endif /* _freeAML_cg_h_ */
