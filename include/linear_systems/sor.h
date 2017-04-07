#ifndef _freeAML_sor_h_
#define _freeAML_sor_h_


#include "generic_iter.h"


namespace aml
{


template < class T >
class sor: public generic_iter< T >
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
	 * @param omega iteration parameter used in the SOR method
	 */
	sor (const size_t max_iter, const T& residual_tol, const T& omega);


	/***********************************************************************
	 *
	 *	SOLVE/ITERATE FUNCTIONS
	 *
	 **********************************************************************/

	/**
	 * @brief solves a linear system using the SOR method
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @return true if the residual tolerance could be achieved within the
	 *         maximum number of iterations allowed, false otherwise
	 */
	template < class _mat, class _vec >
	bool solve (const _mat& A, _vec& x, const _vec& b);


	/***********************************************************************
	 *
	 *	CONFIGURATION FUNCTIONS
	 *
	 **********************************************************************/

	/** @brief gets the value of omega set for the SOR solver */
	T omega () const;

	/** @brief sets the value of omega for the SOR solver */
	void set_omega (const T& omega);


private:

	T	m_omega;

}; /* end of class sor */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
sor< T >::sor(const size_t max_iter,
              const T& residual_tol,
              const T& omega):

        generic_iter< T >(max_iter, residual_tol), m_omega(omega)
{
	FREEAML_ASSERT(omega > (T) 0);

	/* nothing needs to  be done here */
}


template < class T >
template < class _mat, class _vec >
bool sor< T >::solve (const _mat& A, _vec& x, const _vec& b)
{
	FREEAML_ASSERT(generic_iter< T >::check_dimensions(A,x,b));

	_vec x0 = x;

	size_t num_iter = 0;

	while (num_iter < gen_it_type::max_num_iterations())
	{
		for (size_t i = 0; i < x.size(); i++)
		{
			/* if A(i,i) = 0, the SOR method cannot be used */
			if (A(i,i) == (T) 0)
			{
				return false;
			}

			T c = b[i];

			for (size_t j = 0; j < i; j++)
			{
				c -= A(i,j) * x[j];
			}

			for (size_t j = i+1; j < A.num_rows(); j++)
			{
				c -= A(i,j) * x0[j];
			}

			x[i] = ((T)1 - omega())*x0[i] + omega()*(c / A(i,i));
		}

		num_iter++;

		/* if the residual is within the maximum tolerance, stop */
		if ((b - A*x).l2_norm() <= gen_it_type::residual_tolerance())
		{
			gen_it_type::set_num_iterations(num_iter);
			return true;
		}

		x0 = x;
	}

	gen_it_type::set_num_iterations(num_iter);
	return false;
}


template < class T >
T sor< T >::omega () const
{
	return m_omega;
}


template < class T >
void sor< T >::set_omega (const T& omega)
{
	FREEAML_ASSERT(omega > (T) 0);

	m_omega = omega;
}


} /* end of namespace aml */

#endif /* _freeAML_sor_h_ */
