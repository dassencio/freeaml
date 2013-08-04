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


#ifndef _freeAML_plu_pp_rec_h_
#define _freeAML_plu_pp_rec_h_


#include "data_structures/vec.h"
#include "data_structures/sparse_mat.h"


namespace aml
{


template < class T >
class plu_pp_rec
{

public:

	/***********************************************************************
	 *
	 *	FACTORIZE FUNCTIONS
	 *
	 **********************************************************************/

	/**
	 * @brief computes the PLU factorization of an input matrix
	 * @param A the matrix which must be factorized
	 * @param P on success, a permutation matrix
	 * @param L on success, a lower triangular matrix with all diagonal
	 *        elements set to one (1)
	 * @param U on success, an upper triangular nonsingular matrix
	 * @return true if the factorization A = PLU could be built, false otherwise
	 * @note the PLU factorization exists only if A is invertible
	 * @note the recursive procedure is slow (use other PLU factorization
	 *       methods if you need faster performance)
	 * @note partial pivoting is used when building the PLU factorization
	 */
	template < class _mat1, class _mat2,
	           class _mat3, class _mat4 >
	bool factorize (_mat1  A,
	                _mat2& P,
	                _mat3& L,
	                _mat4& U) const;

}; /* end of class plu_pp_rec */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
template < class _mat1, class _mat2,
           class _mat3, class _mat4 >
bool plu_pp_rec< T >::factorize (_mat1  A,
                                 _mat2& P,
                                 _mat3& L,
                                 _mat4& U) const
{
	FREEAML_ASSERT(A.is_square());

	size_t n = A.num_rows();

	/* assign the correct dimensions to P, L and U */
	P.resize(n,n);
	L.resize(n,n);
	U.resize(n,n);

	P.zero_fill();
	L.zero_fill();
	U.zero_fill();

	/* if A is a 1×1 matrix */
	if (n == 1)
	{
		P(0,0) = (T) 1;
		L(0,0) = (T) 1;
		U(0,0) = A(0,0);

		/* if the only entry of A is zero, A is not invertible */
		if (A(0,0) == (T) 0)
		{
			return false;
		}

		return true;
	}

	/* find the largest (in magnitude) element on the first column of A */
	size_t p = 0;
	size_t q = 0;

	while (q < n)
	{
		if (std::abs(A(p,0)) < std::abs(A(q,0)))
		{
			p = q;
		}
		q++;
	}

	/* if the first column of A is  A is not invertible */
	if (A(p,0) == (T) 0)
	{
		return false;
	}

	/* define the matrix Q which permutes the first and p-th rows of A */
	sparse_mat< T > Q(n,n);

	for (size_t i = 1; i < n; i++)
	{
		if (i != p)
		{
			Q(i,i) = (T) 1;
		}
	}

	Q(0,p) = (T) 1;
	Q(p,0) = (T) 1;

	_mat1 B = Q*A;

	/* initialize the matrices A1,P1,L1,U1 for the next level of the
	 * PLU factorization */
	_mat1 A1(n-1,n-1);
	_mat2 P1(n-1,n-1);
	_mat3 L1(n-1,n-1);
	_mat4 U1(n-1,n-1);

	vec< T > U01(n-1, (T)0);
	vec< T > L10(n-1, (T)0);

	/* define L10 and U01 */
	for (size_t i = 1; i < n; i++)
	{
		U01[i-1] = B(0,i);
		L10[i-1] = B(i,0);
	}

	/* define A1 */
	for (size_t i = 1; i < n; i++)
	{
		for (size_t j = 1; j < n; j++)
		{
			A1(i-1,j-1) = B(i,j) - L10[i-1]*U01[j-1] / B(0,0);
		}
	}

	/* continue the PLU factorization recursively */
	if (!factorize(A1, P1, L1, U1))
	{
		return false;
	}

	_mat1 R(n,n);

	R.zero_fill();

	R(0,0) = (T) 1;
	L(0,0) = (T) 1;
	U(0,0) = B(0,0);

	vec< T > P1tL10 = P1.transpose() * L10;

	/* construct the matrices U,L,R at this level of the PLU factorization */
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			if (i > 0 && j > 0)
			{
				R(i,j) = P1(i-1,j-1);
				L(i,j) = L1(i-1,j-1);
				U(i,j) = U1(i-1,j-1);
			}
			else if (i > 0 && j == 0)
			{
				L(i,0) = P1tL10[i-1] / B(0,0);
			}
			else if (j > 0 && i == 0)
			{
				U(0,j) = U01[j-1];
			}
		}
	}

	/* finish constructing the permutation matrix for this level of the
	 * PLU factorization */
	P = Q.transpose() * R;

	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_plu_pp_rec_h_ */
