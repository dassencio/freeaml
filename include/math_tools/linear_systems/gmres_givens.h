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


#ifndef _freeAML_gmres_givens_h_
#define _freeAML_gmres_givens_h_


#include "general/math.h"

#include "data_structures/vec.h"
#include "data_structures/mat.h"

#include "generic_iter.h"


namespace aml
{


template < class T >
class gmres_givens: public generic_iter< T >
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
	gmres_givens (const size_t max_iter, const T& residual_tol);


	/***********************************************************************
	 *
	 *	SOLVE FUNCTION
	 *
	 **********************************************************************/

	/**
	 * @brief solves a linear system using the GMRES method with Givens
	 *        rotations
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @return true if the residual tolerance could be achieved within the
	 *         maximum number of iterations allowed, false otherwise
	 * @note this is a fast version of GMRES (gmres_givens_slow is slower)
	 */
	template < class _mat, class _vec >
	bool solve (const _mat& A, _vec& x, const _vec& b);

}; /* end of class gmres_givens */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
gmres_givens< T >::gmres_givens(const size_t max_iter,
                                const T& residual_tol):

        gen_it_type(max_iter, residual_tol)
{
	/* nothing needs to  be done here */
}


template < class T >
template < class _mat, class _vec >
bool gmres_givens< T >::solve (const _mat& A, _vec& x, const _vec& b)
{
	FREEAML_ASSERT(gen_it_type::check_dimensions(A,x,b));

	size_t n = A.num_rows();

	mat< T > Q(1, n, (T)0);
	mat< T > H(1, n, (T)0);

	T norm_b = b.l2_norm();

	/* build q_0, the basis of the Krylov subspace K^1 = {b} */
	for (size_t j = 0; j < Q.num_cols(); j++)
	{
		Q(0,j) = b[j] / norm_b;
	}

	/* a constant reference to Q, used for accessing its rows directly */
	const mat< T >& c_Q = Q;

	vec< T > f(n, (T)0);

	/* f = (norm_b, 0, ... 0) */
	f[0] = norm_b;

	size_t num_iter = 0;

	/*
	 * c and s will hold the parameters needed to define the successive
	 * Givens rotations applied (not all elements of these vectors will
	 * be used if convergence happens in less than n iterations)
	 */
	vec< T > c(n, (T)0);
	vec< T > s(n, (T)0);

	while (num_iter < gen_it_type::max_num_iterations())
	{
		size_t k = num_iter;

		/* if more rows for Q and H must be allocated */
		if (std::min(k+2,n) > Q.num_rows())
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

		/* z will be used to compute q_k+1 if we have not yet built
		 * the largest Krylov subspace (the one with n dimensions) */
		_vec z = A * c_Q[k];

		/* project out the components of z parallel to {q_0,q_1,...,q_k}
		 * and build h_k (the k-th column of H) */
		for (size_t i = 0; i <= k; i++)
		{
			H(i,k) = c_Q[i] * z;

			z -= H(i,k) * c_Q[i];
		}

		/* apply the previous Givens rotations that kept H upper
		 * triangular (H here being the theoretical H) */
		for (size_t i = 0; i < k; i++)
		{
			givens_rotation(H(i,k), H(i+1,k), c[i], s[i]);
		}

		/* If one more Givens rotations must be applied */
		if (k+1 < n)
		{
			T norm_z = z.l2_norm();

			for (size_t j = 0; j < Q.num_cols(); j++)
			{
				Q(k+1,j) = z[j] / norm_z;
			}

			H(k+1,k) = norm_z;

			/* compute the parameters for the Givens rotation
			 * matrix G(k,k+1) that zeros out H(k+1,k) */
			givens(H(k,k), H(k+1,k), c[k], s[k]);

			/* keep H (the theoretical H) upper triangular */
			givens_rotation(H(k,k), H(k+1,k), c[k], s[k]);

			/* apply also the Givens rotation to the rhs f */
			givens_rotation(f[k], f[k+1], c[k], s[k]);
		}

		/* if H(k,k) = 0, A is not invertible */
		if (H(k,k) == (T) 0)
		{
			break;
		}

		num_iter++;

		/* if we can no longer iterate or if the residual is within
		 * the maximum tolerance */
		if (k == n-1 || std::abs(f[k+1]) <= gen_it_type::residual_tolerance())
		{
			vec< T > lambda(k+1, (T)0);

			/* solve H_{0:k,0:k}*lambda = f */
			size_t i = k+1;

			while (i-- > 0)
			{
				T y = (T) 0;

				for (size_t j = i+1; j <= k; j++)
				{
					y += H(i,j) * lambda[j];
				}

				lambda[i] = (f[i] - y) / H(i,i);
			}

			/* build the problem solution x = Q*lambda (keep in
			 * mind Q is the theoretical Q^t) */
			x.zero_fill();

			for (size_t i = 0; i <= k; i++)
			{
				x += lambda[i] * c_Q[i];
			}

			if (num_iter < gen_it_type::max_num_iterations() ||
			        (b - A*x).l2_norm() <= gen_it_type::residual_tolerance())
			{
				gen_it_type::set_num_iterations(num_iter);
				return true;
			}
		}
	}

	gen_it_type::set_num_iterations(num_iter);
	return false;
}


} /* end of namespace aml */

#endif /* _freeAML_gmres_givens_h_ */
