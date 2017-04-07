#ifndef _freeAML_bicgstab_h_
#define _freeAML_bicgstab_h_


#include "generic_iter.h"


namespace aml
{


template < class T >
class bicgstab: public generic_iter< T >
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
	bicgstab (const size_t max_iter, const T& residual_tol);


	/***********************************************************************
	 *
	 *	SOLVE FUNCTION
	 *
	 **********************************************************************/

	/**
	 * @brief solves a linear system using the biconjugate gradient
	 *        stabilized method
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @return true if the residual tolerance could be achieved within the
	 *         maximum number of iterations allowed, false otherwise
	 */
	template < class _mat, class _vec >
	bool solve (const _mat& A, _vec& x, const _vec& b);

}; /* end of class bicgstab */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
bicgstab< T >::bicgstab(const size_t max_iter, const T& residual_tol):

        gen_it_type(max_iter, residual_tol)
{
	/* nothing needs to  be done here */
}


template < class T >
template < class _mat, class _vec >
bool bicgstab< T >::solve (const _mat& A, _vec& x, const _vec& b)
{
	FREEAML_ASSERT(gen_it_type::check_dimensions(A,x,b));

	/* r = residual vector */
	_vec r = b - A*x;

	_vec u = r;

	T rho   = (T) 1;
	T alpha = (T) 1;
	T omega = (T) 1;

	_vec p(b.size());
	_vec v(b.size());

	p.zero_fill();
	v.zero_fill();

	size_t num_iter = 0;

	while (num_iter <= gen_it_type::max_num_iterations())
	{
		T old_rho = rho;

		rho = u*r;

		T beta = (rho * alpha) / (old_rho * omega);

		p = r + beta*(p - omega*v);
		v = A*p;

		alpha = rho / (u*v);

		_vec s = r - alpha*v;
		_vec t = A*s;

		omega = (t*s) / (t*t);

		x += alpha*p + omega*s;

		r = s - omega*t;

		num_iter++;

		/* if the residual is within the maximum tolerance, stop */
		if (r.l2_norm() <= gen_it_type::residual_tolerance())
		{
			gen_it_type::set_num_iterations(num_iter);
			return true;
		}
	}

	gen_it_type::set_num_iterations(num_iter);
	return false;
}


} /* end of namespace aml */

#endif /* _freeAML_bicgstab_h_ */
