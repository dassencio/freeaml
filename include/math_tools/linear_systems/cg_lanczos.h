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


#ifndef _freeAML_cg_lanczos_h_
#define _freeAML_cg_lanczos_h_


#include "generic_iter.h"


namespace aml
{


template < class T >
class cg_lanczos: public generic_iter< T >
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
	cg_lanczos (const size_t max_iter, const T& residual_tol);


	/***********************************************************************
	 *
	 *	SOLVE FUNCTION
	 *
	 **********************************************************************/

	/**
	 * @brief solves a linear system using the Lanczos version of the
	 *        conjugate gradient method
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @return true if the residual tolerance could be achieved within the
	 *         maximum number of iterations allowed, false otherwise
	 */
	template < class _mat, class _vec >
	bool solve (const _mat& A, _vec& x, const _vec& b);

}; /* end of class cg_lanczos */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
cg_lanczos< T >::cg_lanczos(const size_t max_iter, const T& residual_tol):

        gen_it_type(max_iter, residual_tol)
{
	/* nothing needs to  be done here */
}


template < class T >
template < class _mat, class _vec >
bool cg_lanczos< T >::solve (const _mat& A, _vec& x, const _vec& b)
{
	FREEAML_ASSERT(gen_it_type::check_dimensions(A,x,b));

	/* in this algorithm, the initial guess is the zero vector */
	x.zero_fill();

	_vec v = b;

	/* \beta_{k-1} and \beta_k, k is the iteration number */
	T beta1 = (T) 0;
	T beta2 = v.l2_norm();

	/* d_{k-1} and d_k */
	T delta1 = (T) 0;
	T delta2 = (T) 0;

	/* vectors q_{k-1} and q_k */
	_vec q1(x.size());
	_vec q2(x.size());

	q1.zero_fill();
	q2.zero_fill();

	_vec c(x.size(), (T) 0);

	/* \rho_k */
	T rho = (T) 0;

	size_t num_iter = 0;

	while (num_iter < gen_it_type::max_num_iterations())
	{
		/* compute the Lanczos vector q_k, set the ground for q_{k+1} */

		q1 = q2;
		q2 = v / beta2;
		v  = A*q2;

		T alpha = q2*v;

		v -= alpha*q2 + beta2*q1;

		beta1 = beta2;
		beta2 = v.l2_norm();

		/* LDL^t decomposition with forward substitution */
		if (num_iter == 0)
		{
			delta2 = alpha;

			c = q2;

			rho = beta1 / alpha;

			x = rho*q2;
		}
		else
		{
			delta1 = delta2;

			T mu = beta1 / delta1;

			delta2 = alpha - mu*beta1;

			c = q2 - mu*c;

			rho = -(mu * delta1 / delta2) * rho;

			x += rho*c;
		}

		num_iter++;

		/* if the residual is within the maximum tolerance, stop */
		if ((b - A*x).l2_norm() <= gen_it_type::residual_tolerance())
		{
			gen_it_type::set_num_iterations(num_iter);
			return true;
		}
	}

	gen_it_type::set_num_iterations(num_iter);
	return false;
}


} /* end of namespace aml */

#endif /* _freeAML_cg_lanczos_h_ */
