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


#ifndef _freeAML_qr_house_h_
#define _freeAML_qr_house_h_


#include <cmath>

#include <general/math.h>

#include <data_structures/vec.h>
#include <data_structures/mat.h>


namespace aml
{


template < class T >
class qr_house
{

public:

	/**
	 * @brief computes the QR factorization of a matrix using Householder
	 *        transformations
	 * @param A the matrix which must be factorized
	 * @param Q on success, the matrix Q of the QR factorization
	 * @param R on success, the matrix R of the QR factorization
	 * @return true (this function always succeeds)
	 * @note at the end Q will be orthogonal and R will be upper
	 *       triangular with nonnegative diagonal values
	 */
	template < class _mat1, class _mat2, class _mat3 >
	bool factorize (_mat1 A, _mat2& Q, _mat3& R) const;

}; /* end of class qr_house */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
template < class _mat1, class _mat2, class _mat3 >
bool qr_house< T >::factorize (_mat1 A, _mat2& Q, _mat3& R) const
{
	size_t m = A.num_rows();
	size_t n = A.num_cols();

	/* P will be the product of all Householder transformations applied */
	mat< T > P = dense_identity_matrix< T >(m);

	size_t d = std::min(m,n);

	/* assign the correct dimensions to Q and R */
	Q.resize(m, d);
	R.resize(d, n);

	Q.zero_fill();
	R.zero_fill();

	vec< T > w(n, (T) 0);
	vec< T > z(m, (T) 0);

	/* for each column p of A which must have its subdiagonal
	 * elements zeroed out */
	for (size_t p = 0; p < d; p++)
	{
		size_t M = m-p;
		size_t N = n-p;

		/* define u := A(p:m-1,p) */
		vec< T > u(M, (T) 0);

		for (size_t i = 0; i < M; i++)
		{
			u[i] = A(i+p,p);
		}

		/* define the Householder transformation for u which
		 * projects it onto e_0 = (1,0,...0) */
		vec< T > v;

		T b,c = (T) 0;

		householder(u,v,b,c,0);

		/* If b = 0, no Householder transformation is necessary */
		if (b != (T) 0)
		{
			/* build w = (v^t)A(p:m-1,p:n-1) and
			 * z = (v^t)P(p:m-1,0:m-1) */

			w.zero_fill();
			z.zero_fill();

			for (size_t k = 0; k < M; k++)
			{
				for (size_t j = 0; j < N; j++)
				{
					w[j] += v[k] * A(k+p,j+p);
				}

				for (size_t j = 0; j < m; j++)
				{
					z[j] += v[k] * P(k+p,j);
				}
			}

			v *= (-b);

			/* A(p:m-1,p:n-1) <--- -c*(I - bvv^t)A(p:m-1,p:n-1) */
			/* P(p:m-1,0:m-1) <--- -c*(I - bvv^t)P(p:m-1,0:m-1) */
			for (size_t i = 0; i < M; i++)
			{
				/* update the rows p to n-1 of A */
				for (size_t j = 0; j < N; j++)
				{
					A(i+p,j+p) = -c*(A(i+p,j+p) + v[i]*w[j]);
				}

				/* update the rows p to m-1 of P */
				for (size_t j = 0; j < m; j++)
				{
					P(i+p,j) = -c*(P(i+p,j) + v[i]*z[j]);
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

#endif /* _freeAML_qr_house_h_ */
