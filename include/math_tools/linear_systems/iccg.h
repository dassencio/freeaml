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


#ifndef _freeAML_iccg_h_
#define _freeAML_iccg_h_


#include <cmath>

#include "generic_iter.h"


namespace aml
{


template < class T >
class iccg: public generic_iter< T >
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
	iccg (const size_t max_iter, const T& residual_tol);


	/***********************************************************************
	 *
	 *	SOLVE FUNCTION
	 *
	 **********************************************************************/

	/**
	 * @brief solves a linear system using the conjugate gradient method
	 *        with the incomplete Cholesky preconditioner
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @return true if the residual tolerance could be achieved within the
	 *         maximum number of iterations allowed, false otherwise
	 */
	template < class _mat, class _vec >
	bool solve (const _mat& A, _vec& x, const _vec& b);


private:

	/**
	 * @brief builds the lower triangular matrix K which approximates the
	 *        Cholesky factorization of an input matrix (KK^t)
	 * @param A the input matrix
	 * @param on success, the lower triangular matrix K such that A is
	 *        "approximately equal to" KK^t
	 */
	template < class _mat >
	bool build_K (const _mat& A, _mat& K) const;

	/**
	 * @brief solves the ICCG preconditioner equation KK^t*z = r
	 * @param K the lower triangular matrix which approximates the Cholesky
	 *        factorization of the system matrix (A)
	 * @param z a vector on which the solution to the preconditioner
	 *        equation will be written (if found)
	 * @param r the right-hand side vector (in our case, a residual vector)
	 */
	template < class _mat, class _vec >
	bool solve_preconditioner_equation (const _mat& K,
	                                          _vec& z,
	                                    const _vec& r) const;

}; /* end of class iccg */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
iccg< T >::iccg(const size_t max_iter, const T& residual_tol):

        gen_it_type(max_iter, residual_tol)
{
	/* nothing needs to  be done here */
}


template < class T >
template < class _mat, class _vec >
bool iccg< T >::solve (const _mat& A, _vec& x, const _vec& b)
{
	FREEAML_ASSERT(gen_it_type::check_dimensions(A,x,b));

	size_t n = A.num_rows();

	_mat K(n,n);

	K.set_all_values_to_zero();

	/* the incomplete Cholesky preconditioner is KK^t */
	if (!build_K(A,K))
	{
		gen_it_type::set_num_iterations(0);
		return false;
	}

	/* r = residual vector */
	_vec r = b - A*x;

	_vec z(n);

	/* do z = (KK^t)^(-1)r */
	if (!solve_preconditioner_equation(K,z,r))
	{
		gen_it_type::set_num_iterations(0);
		return false;
	}

	_vec p = z;
	_vec q = A*z;

	size_t num_iter = 0;

	while (num_iter <= gen_it_type::max_num_iterations())
	{
		T gamma = r*z;
		T alpha = gamma / (p*q);

		x += alpha * p;
		r -= alpha * q;

		/* do z = (KK^t)^(-1)r */
		if (!solve_preconditioner_equation(K,z,r))
		{
			break;
		}

		T beta = (r*z) / gamma;

		p = z + (beta * p);
		q = A*z + (beta * q);

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


template < class T >
template < class _mat >
bool iccg< T >::build_K (const _mat&  A, _mat& K) const
{
	typedef typename _mat::const_iterator const_iterator;

	size_t n = A.num_rows();

	/* const reference to K (to read it without changing it) */
	const _mat& c_K = K;

	for (size_t j = 0; j < n; j++)
	{
		if (A(j,j) != (T) 0)
		{
			K(j,j) = A(j,j);

			for (const_iterator k = c_K[j].begin(); k->first < j; k++)
			{
				K(j,j) -= k->second * k->second;
			}
		}

		K(j,j) = std::sqrt(K(j,j));

		/* if K(j,j) = 0, then A is not invertible */
		if (K(j,j) == (T) 0)
		{
			return false;
		}

		for (size_t i = j+1; i < n; i++)
		{
			if (A(i,j) != (T) 0)
			{
				K(i,j) = A(i,j);

				for (const_iterator k = c_K[i].begin(); k->first < j; k++)
				{
					K(i,j) -= k->second * c_K(j,k->first);
				}

				K(i,j) /= K(j,j);
			}
		}
	}

	return true;
}


template < class T >
template < class _mat, class _vec >
bool iccg< T >::solve_preconditioner_equation (const _mat& K,
                                                     _vec& z,
                                               const _vec& r) const
{
	typedef typename _mat::const_iterator const_iterator;

	/* K is an n × n matrix */
	size_t n = K.num_rows();

	z.set_all_values_to_zero();

	_vec w(n);

	w.set_all_values_to_zero();

	const _mat Kt = K.transpose();

	/* perform forward substitution to solve Kw = b, with w = (K^t)z */
	for (size_t i = 0; i < n; i++)
	{
		T y = (T) 0;

		for (const_iterator j = K[i].begin(); j->first < i; j++)
		{
			y += j->second * w[j->first];
		}

		w[i] = (r[i] - y) / K(i,i);
	}

	/* perform backward substitution to solve (K^t)z = w */
	size_t i = n;

	while (i-- > 0)
	{
		T y = (T) 0;

		for (const_iterator j = Kt[i].begin(); j != Kt[i].end(); j++)
		{
			y += j->second * z[j->first];
		}

		z[i] = (w[i] - y) / Kt(i,i);
	}

	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_iccg_h_ */
