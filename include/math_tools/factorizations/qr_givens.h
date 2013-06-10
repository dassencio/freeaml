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


#ifndef _freeAML_qr_givens_h_
#define _freeAML_qr_givens_h_


#include <cmath>

#include "general/math.h"

#include "data_structures/mat.h"


namespace aml
{


template < class T >
class qr_givens
{

public:

	/**
	 * @brief computes the QR factorization of a matrix using Givens
	 *        rotations
	 * @param A the matrix which must be factorized
	 * @param Q on success, the matrix Q of the QR factorization
	 * @param R on success, the matrix R of the QR factorization
	 * @return true (this function always succeeds)
	 * @note at the end Q will be orthogonal and R will be upper
	 *       triangular with nonnegative diagonal values
	 */
	template < class _mat1, class _mat2, class _mat3 >
	bool factorize (_mat1 A, _mat2& Q, _mat3& R) const;

}; /* end of class qr_givens */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
template < class _mat1, class _mat2, class _mat3 >
bool qr_givens< T >::factorize (_mat1 A, _mat2& Q, _mat3& R) const
{
	size_t m = A.num_rows();
	size_t n = A.num_cols();

	/* P will be the product of all Givens rotations applied */
	mat< T > P = identity_matrix< T >(m);

	size_t d = std::min(m,n);

	/* assign the correct dimensions to Q and R */
	Q.resize(m, d);
	R.resize(d, n);

	Q.set_all_values_to_zero();
	R.set_all_values_to_zero();

	/* for each column p of A which must have its subdiagonal
	 * elements zeroed out */
	for (size_t p = 0; p < d; p++)
	{
		/* for each element on the p-th column of A which must be
		 * zeroed out */
		for (size_t i = p+1; i < m; i++)
		{
			/* if A(i,p) is not already zero */
			if (A(i,p) != (T) 0)
			{
				/* define the Givens rotation matrix parameters
				 * c and s for zeroing out A(i,p) */
				T c = (T) 0;
				T s = (T) 0;

				givens(A(p,p), A(i,p), c, s);

				/* apply the Givens rotation to A */
				for (size_t j = p; j < n; j++)
				{
					T a = A(p,j);
					T b = A(i,j);

					A(p,j) = c*a - s*b;
					A(i,j) = s*a + c*b;
				}

				/* apply the Givens rotation to P */
				for (size_t j = 0; j < m; j++)
				{
					T a = P(p,j);
					T b = P(i,j);

					P(p,j) = c*a - s*b;
					P(i,j) = s*a + c*b;
				}
			}
		}
	}

	/* build Q from the first m rows of P^t */
	for (size_t i = 0; i < m; i++)
	{
		for (size_t j = 0; j < d; j++)
		{
			Q(i,j) = P(j,i);
		}
	}

	/* at this point A is upper triangular, so build R from its
	 * first d rows */
	for (size_t i = 0; i < d; i++)
	{
		for (size_t j = i; j < n; j++)
		{
			R(i,j) = A(i,j);
		}
	}

	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_qr_givens_h_ */