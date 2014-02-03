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


#ifndef _freeAML_gmres_house_h_
#define _freeAML_gmres_house_h_


#include <data_structures/vec.h>
#include <data_structures/mat.h>

#include <factorizations/qr_house.h>
#include <least_squares/qr.h>

#include "generic_iter.h"


namespace aml
{


template < class T >
class gmres_house: public generic_iter< T >
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
	gmres_house (const size_t max_iter, const T& residual_tol);


	/***********************************************************************
	 *
	 *	SOLVE FUNCTION
	 *
	 **********************************************************************/

	/**
	 * @brief solves a linear system using the GMRES method with Householder
	 *        transformations
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
	 * @brief performs the k-th step of the Arnoldi algorithm for a matrix
	 * @param A the system matrix (will be decomposed)
	 * @param H on success, the upper Hessenberg matrix associated with the
	 *        Hessenberg decomposition of A on the k-th step of the Arnoldi
	 *        algorithm
	 * @param Q on success, the transpose of the orthogonal matrix associated
	 *        with the Heissenberg decomposition of A on the k-th step of
	 *        the Arnoldi algorithm
	 * @param b the right-hand side vector
	 * @param k the Arnoldi algorithm step number
	 * @return true if the step succeeded, false otherwise
	 */
	template < class _mat1, class _mat2, class _vec >
	bool arnoldi_iteration (const _mat1&  A,
	                              _mat2&  H,
	                              _mat2&  Q,
	                        const _vec&   b,
	                        const size_t  k) const;

}; /* end of class gmres_house */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
gmres_house< T >::gmres_house(const size_t max_iter, const T& residual_tol):

        gen_it_type(max_iter, residual_tol)
{
	/* nothing needs to  be done here */
}


template < class T >
template < class _mat, class _vec >
bool gmres_house< T >::solve (const _mat& A, _vec& x, const _vec& b)
{
	FREEAML_ASSERT(gen_it_type::check_dimensions(A,x,b));

	size_t n = A.num_rows();

	/*
	 * Q and H are the matrices which define the Hessenberg decomposition
	 * of A; Q is actually the transpose of the theoretical Q (we allocate
	 * space for these matrices as necessary because they are dense)
	 */
	mat< T > Q(1, n, (T)0);
	mat< T > H(1, n, (T)0);

	/* a constant reference to Q, used for accessing its rows directly */
	const mat< T >& c_Q = Q;

	/* QR factorization class with Householder transformations */
	qr_house< T > fac;

	/* general least squares solver with QR factorization */
	qr< T > lss;

	T norm_b = b.l2_norm();

	size_t num_iter = 0;

	while (num_iter < gen_it_type::max_num_iterations())
	{
		size_t k = num_iter;

		/*
		 * s defines the number of rows of Hk and fk (see below); it
		 * is the number of basis vectors q_i built at the end this
		 * iteration
		 */
		size_t s = std::min(k+2, n);

		/* if more rows for Q and H must be allocated */
		if (s > Q.num_rows())
		{
			/*
			 * in order to avoid wasting CPU cycles, we double the
			 * number of rows allocated for Q and H whenever
			 * necessary instead of allocating one row at a time
			 */
			size_t m = std::min(n, 2*Q.num_rows());

			Q.resize(m, n, (T)0);
			H.resize(m, n, (T)0);
		}

		/* execute the next step of the Arnoldi algorithm */
		if (!arnoldi_iteration(A, H, Q, b, k))
		{
			break;
		}

		vec< T > y(k+1, (T) 0);

		/* Hk(0:s-1, 0:k) <--- H(0:s-1, 0:k) */
		mat< T > Hk = H.submatrix(0, s-1, 0, k);

		vec< T > fk(s, (T) 0);

		/* fk = (norm_b, 0, ... 0) */
		fk[0] = norm_b;

		/* find y which minimizes |Hk*y - fk|_2; y is in
		 * K^{k+1} = span{b,Ab,...,A^kb} */
		lss.solve(Hk, y, fk, fac);

		/* update the approximate solution x */
		x.zero_fill();

		for (size_t i = 0; i <= k; i++)
		{
			x += y[i] * c_Q[i];
		}

		num_iter++;

		/* if the residual is within the maximum tolerance, stop */
		if ((b - A*x).linf_norm() <= gen_it_type::residual_tolerance())
		{
			gen_it_type::set_num_iterations(num_iter);
			return true;
		}
	}

	gen_it_type::set_num_iterations(num_iter);
	return false;
}


template < class T >
template < class _mat1, class _mat2, class _vec >
bool gmres_house< T >::arnoldi_iteration (const _mat1&  A,
                                                _mat2&  H,
                                                _mat2&  Q,
                                          const _vec&   b,
                                          const size_t  k) const
{
	size_t n = A.num_rows();

	if (k == 0)
	{
		T norm_b = b.l2_norm();

		if (norm_b > (T) 0)
		{
			for (size_t j = 0; j < b.size(); j++)
			{
				Q(0,j) = b[j] / norm_b;
			}
		}
		else
		{
			return false;
		}
	}

	const mat< T >& c_Q = Q;

	/* z will later be used to compute q_k+1, the (k+1)-th row of Q */
	vec< T > z = A * c_Q[k];

	/* project out q_0,...,q_k from z to get an orthogonal basis of
	 * span{b,Ab,...A^k*b} */
	for (size_t i = 0; i <= k; i++)
	{
		H(i,k) = c_Q[i] * z;

		z -= H(i,k) * c_Q[i];
	}

	/* finish building q_k+1 and get H_{k+1,k} */
	if (k < n-1)
	{
		H(k+1,k) = z.l2_norm();

		if (H(k+1,k) == (T) 0)
		{
			return false;
		}

		/* normalize q_k+1 */
		for (size_t j = 0; j < z.size(); j++)
		{
			Q(k+1,j) = z[j] / H(k+1,k);
		}
	}

	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_gmres_house_h_ */
