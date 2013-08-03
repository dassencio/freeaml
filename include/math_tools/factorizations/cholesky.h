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


#ifndef _freeAML_cholesky_h_
#define _freeAML_cholesky_h_


#include <cmath>

#include "general/debug.h"


namespace aml
{


template < class T >
class cholesky
{

public:

	/**
	 * @brief computes the Cholesky factorization of an input matrix
	 * @param A the matrix which must be factorized
	 * @param L on success, a lower triangular matrix with positive
	 *        diagonal entries such that A = LL^t is the Cholesky
	 *        factorization of A
	 * @return true if the factorization A = LL^t could be built, false
	 *         otherwise
	 */
	template < class _mat1, class _mat2 >
	bool factorize (const _mat1&  A, _mat2& L) const;

}; /* end of class cholesky */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
template < class _mat1, class _mat2 >
bool cholesky< T >::factorize (const _mat1&  A, _mat2& L) const
{
	/* A must be real symmetric positive definite */
	FREEAML_ASSERT(A.is_symmetric());

	size_t n = A.num_rows();

	/* assign the correct dimensions to L */
	L.resize(n,n);

	L.zero_fill();

	for (size_t j = 0; j < n; j++)
	{
		L(j,j) = A(j,j);

		for (size_t k = 0; k < j; k++)
		{
			L(j,j) -= L(j,k) * L(j,k);
		}

		/* if L(j,j) <= 0, A is not positive definite and therefore
		 * has no Cholesky factorization */
		if (L(j,j) > (T) 0)
		{
			L(j,j) = std::sqrt(L(j,j));
		}
		else
		{
			return false;
		}

		for (size_t i = j+1; i < n; i++)
		{
			L(i,j) = A(i,j);

			for (size_t k = 0; k < j; k++)
			{
				L(i,j) -= L(i,k) * L(j,k);
			}

			L(i,j) /= L(j,j);
		}
	}
	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_cholesky_h_ */
