/*******************************************************************************
 *
 * Copyright (c) 2013, Diego Assêncio
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ******************************************************************************/


#ifndef _freeAML_weighted_jacobi_h_
#define _freeAML_weighted_jacobi_h_


#include "generic_iter.h"


namespace aml
{


template < class T >
class weighted_jacobi: public generic_iter< T >
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
	 * @param omega iteration parameter used in the weighted Jacobi
	 *        method (default: 2/3)
	 */
	weighted_jacobi (const size_t max_iter,
	                 const T& residual_tol,
	                 const T& omega = (T)2/(T)3);


	/***********************************************************************
	 *
	 *	SOLVE/ITERATE FUNCTIONS
	 *
	 **********************************************************************/

	/**
	 * @brief solves a linear system using the weighted Jacobi method
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @return true if the residual tolerance could be achieved within the
	 *         maximum number of iterations allowed, false otherwise
	 */
	template < class _mat, class _vec >
	bool solve (const _mat& A, _vec& x, const _vec& b);

	/**
	 * @brief performes n weighted Jacobi iterations for a linear system
	 * @param A the system matrix
	 * @param x the solution vector
	 * @param b the right-hand side vector
	 * @note this function blindly performs n iterations of the weighted
	 *       Jacobi method without taking the residual tolerance into
	 *       consideration
	 */
	template < class _mat, class _vec >
	void iterate (const _mat& A, _vec& x, const _vec& b, const size_t n);


	/***********************************************************************
	 *
	 *	CONFIGURATION FUNCTIONS
	 *
	 **********************************************************************/

	/** @brief gets the value of omega set for the Jacobi solver */
	T omega () const;

	/** @brief sets the value of omega for the Jacobi solver */
	void set_omega (const T& omega);


private:

	T	m_omega;

}; /* end of class weighted_jacobi */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
weighted_jacobi< T >::weighted_jacobi(const size_t max_iter,
                                      const T& residual_tol,
		                      const T& omega):

        generic_iter< T >(max_iter, residual_tol), m_omega(omega)
{
	FREEAML_ASSERT(omega > (T) 0);

	/* nothing needs to  be done here */
}


template < class T >
template < class _mat, class _vec >
bool weighted_jacobi< T >::solve (const _mat& A, _vec& x, const _vec& b)
{
	FREEAML_ASSERT(generic_iter< T >::check_dimensions(A,x,b));

	/* r = residual vector */
	_vec r = b - A*x;

	size_t num_iter = 0;

	while (num_iter < gen_it_type::max_num_iterations())
	{
		for (size_t i = 0; i < x.size(); i++)
		{
			FREEAML_ASSERT(A(i,i) != (T) 0);

			x[i] += omega() * r[i] / A(i,i);
		}

		num_iter++;

		r = b - A*x;

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


template < class T >
template < class _mat, class _vec >
void weighted_jacobi< T >::iterate (const _mat& A,
                                          _vec& x,
                                    const _vec& b,
                                    const size_t n)
{
	FREEAML_ASSERT(n > 0);
	FREEAML_ASSERT(generic_iter< T >::check_dimensions(A,x,b));

	/* r = residual vector */
	_vec r = b - A*x;

	size_t num_iter = 0;

	while (num_iter < n)
	{
		for (size_t i = 0; i < x.size(); i++)
		{
			x[i] += omega() * r[i] / A(i,i);
		}

		num_iter++;

		r = b - A*x;
	}

	generic_iter< T >::set_num_iterations(num_iter);
}


template < class T >
T weighted_jacobi< T >::omega () const
{
	return m_omega;
}


template < class T >
void weighted_jacobi< T >::set_omega (const T& omega)
{
	FREEAML_ASSERT(omega > (T) 0);

	m_omega = omega;
}


} /* end of namespace aml */

#endif /* _freeAML_weighted_jacobi_h_ */