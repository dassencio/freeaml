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


#ifndef _freeAML_bidiagonal_h_
#define _freeAML_bidiagonal_h_


#include "data_structures/vec.h"
#include "data_structures/mat.h"


namespace aml
{


template < class T >
class bidiagonal
{

public:

	/***********************************************************************
	 *
	 *	FACTORIZE FUNCTIONS
	 *
	 **********************************************************************/

	/**
	 * @brief computes the bidiagonal factorization of an input matrix
	 * @param A the matrix which must be factorized
	 * @param U on success, an orthogonal matrix
	 * @param V on success, an orthogonal matrix
	 * @param Z on success, an upper bidiagonal matrix
	 * @return true if the bidiagonal factorization A = UZV^t could be
	 *         built, false otherwise
	 * @note Householder transformations are used to compute U,V,Z
	 */
	template < class _mat1, class _mat2,
	           class _mat3, class _mat4 >
	bool factorize (_mat1  A,
	                _mat2& U,
	                _mat3& V,
	                _mat4& Z) const;

}; /* end of class bidiagonal */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
template < class _mat1, class _mat2,
           class _mat3, class _mat4 >
bool bidiagonal< T >::factorize (_mat1  A,
                                 _mat2& U,
                                 _mat3& V,
                                 _mat4& Z) const
{
	size_t m = A.num_rows();
	size_t n = A.num_cols();

	/* assign the correct dimensions to U and V (U and V will be products
	 * of Householder transformations) */
	U.resize(m,m);
	V.resize(n,n);

	U.zero_fill();
	V.zero_fill();

	/* initialize U to the m×m identity matrix */
	for (size_t i = 0; i < m; i++)
	{
		U(i,i) = (T) 1;
	}

	/* initialize V to the n×n identity matrix */
	for (size_t i = 0; i < n; i++)
	{
		V(i,i) = (T) 1;
	}

	size_t d = std::min(m,n);

	/* for each column p of B which must have the portion below
	 * its subdiagonal zeroed out */
	for (size_t p = 0; p < d; p++)
	{
		size_t M = m-p;
		size_t N = n-p;

		/* let uL := A(p:m-1,p) */
		vec< T > uL(M, (T) 0);

		for (size_t i = 0; i < M; i++)
		{
			uL[i] = A(i+p,p);
		}

		/* define the Householder transformation for uL which
		 * projects it onto e_0 = (1,0,...0) */
		vec< T > vL;

		T b,c = 0;

		householder(uL,vL,b,c,0);

		/* If b = 0, no Householder transformation is necessary */
		if (b != (T) 0)
		{
			mat< T > B(M, N, (T) 0);
			mat< T > C(M, m, (T) 0);

			for (size_t i = 0; i < M; i++)
			{
				/* B <-- (I - bvv^t)A(p:m-1,p:n-1), where v = vL */
				for (size_t j = 0; j < N; j++)
				{
					B(i,j) = A(i+p,j+p);

					for (size_t k = 0; k < M; k++)
					{
						B(i,j) -= b * vL[i] * vL[k] * A(k+p,j+p);
					}
				}

				/* C <-- (I - bvv^t)U(p:m-1,0:m-1), where v = vL */
				for (size_t j = 0; j < m; j++)
				{
					C(i,j) = U(i+p,j);

					for (size_t k = 0; k < M; k++)
					{
						C(i,j) -= b * vL[i] * vL[k] * U(k+p,j);
					}
				}
			}

			for (size_t i = 0; i < M; i++)
			{
				/* A(p:m-1,p:n-1) <--- -c*B */
				for (size_t j = 0; j < N; j++)
				{
					A(i+p,j+p) = -c * B(i,j);
				}

				/* U(p:m-1,0:m-1) <--- -c*C */
				for (size_t j = 0; j < m; j++)
				{
					U(i+p,j) = -c * C(i,j);
				}
			}
		}

		/* if V needs to be updated (the "last" V is just the identity) */
		if (N-1 > 0)
		{
			/* let uR := A(p,p+1:n-1) */
			vec< T > uR(N-1, (T) 0);

			for (size_t j = 0; j < N-1; j++)
			{
				uR[j] = A(p,p+j+1);
			}

			/* define the Householder transformation for uR which
			 * projects it onto e_0 = (1,0,...0) */
			vec< T > vR;

			c = (T) 0;
			b = (T) 0;

			householder(uR,vR,b,c,0);

			/* If b = 0, no Householder transformation is necessary */
			if (b != (T) 0)
			{
				mat< T > D(M, N-1, (T) 0);
				mat< T > E(n, N-1, (T) 0);

				for (size_t j = 0; j < N-1; j++)
				{
					/* D <-- A(p:m-1,p+1:n-1)(I - b(vv^t)^t), where v = vR */
					for (size_t i = 0; i < M; i++)
					{
						D(i,j) = A(p+i,p+j+1);

						for (size_t k = 0; k < N-1; k++)
						{
							D(i,j) -= A(p+i,p+k+1) * b * vR[k] * vR[j];
						}
					}

					/* E = V(0:n-1,p+1:n-1)(I - b(vv^t)^t), where v = vR */
					for (size_t i = 0; i < n; i++)
					{
						E(i,j) = V(i,p+j+1);

						for (size_t k = 0; k < N-1; k++)
						{
							E(i,j) -= V(i,p+k+1) * b * vR[k] * vR[j];
						}
					}
				}

				for (size_t j = 0; j < N-1; j++)
				{
					/* A(p:m-1,p+1:n-1) <--- -c*D */
					for (size_t i = 0; i < M; i++)
					{
						A(p+i,p+j+1) = -c * D(i,j);
					}

					/* V(0:n-1,p+1:n-1) <--- -c*E */
					for (size_t i = 0; i < n; i++)
					{
						V(i,p+j+1) = -c * E(i,j);
					}
				}
			}
		}
	}

	/* Z <-- A (at this point, A is the biadiagonal form of the original A */
	Z = A;

	/*U <--- U^t */
	U = U.transpose();

	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_bidiagonal_h_ */
