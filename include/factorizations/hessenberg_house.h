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


#ifndef _freeAML_hessenberg_house_h_
#define _freeAML_hessenberg_house_h_


#include <data_structures/vec.h>
#include <data_structures/mat.h>


namespace aml
{


template < class T >
class hessenberg_house
{

public:

	/***********************************************************************
	 *
	 *	FACTORIZE FUNCTIONS
	 *
	 **********************************************************************/
	/**
	 * @brief computes the Hessenberg decomposition of an input matrix
	 * @param A the input matrix which must be factorized
	 * @param Q at the end, an orthogonal matrix
	 * @param H at the end, an upper Hessenberg matrix with nonnegative
	 *        values on its first subdiagonal
	 * @return true (this function always succeeds)
	 * @note the Hessenberg decomposition is given by A = QHQ^t; in most
	 *       cases H is a dense matrix
	 */
	template < class _mat1, class _mat2, class _mat3 >
	bool factorize (const _mat1& A, _mat2& Q, _mat3& H) const;

	/**
	 * @brief computes the Hessenberg decomposition of an input matrix
	 * @param A the input matrix which must be factorized
	 * @param H at the end, an upper Hessenberg matrix with nonnegative
	 *        values on its first subdiagonal
	 * @return true (this function always succeeds)
	 * @note the Hessenberg decomposition is given by A = QHQ^t, but this
	 *       function computes only H (Q is not computed); in most cases
	 *       H is a dense matrix
	 */
	template < class _mat1, class _mat2 >
	bool factorize (const _mat1& A, _mat2& H) const;

}; /* end of class hessenberg_house */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
template < class _mat1, class _mat2, class _mat3 >
bool hessenberg_house< T >::factorize (const _mat1& A, _mat2& Q, _mat3& H) const
{
	FREEAML_ASSERT(A.is_square());

	size_t m = A.num_rows();

	/* P will be the product of Householder transformations */
	mat< T > P = dense_identity_matrix< T >(m);

	/* assign the correct dimensions to Q and H */
	Q.resize(m,m);
	H.resize(m,m);

	/* H <-- A */
	for (size_t i = 0; i < A.num_rows(); i++)
	{
		for (size_t j = 0; j < A.num_cols(); j++)
		{
			H(i,j) = A(i,j);
		}
	}

	/* for each column p of H which must have the portion below
	 * its subdiagonal zeroed out */
	for (size_t p = 0; p < m-1; p++)
	{
		size_t M = m-p-1;

		/* let u := H(p+1:m-1,p) */
		vec< T > u(M, (T) 0);

		for (size_t i = 0; i < M; i++)
		{
			u[i] = H(i+p+1,p);
		}

		/* define the Householder transformation for u which
		 * projects it onto e_0 = (1,0,...0) */
		vec< T > v;

		T b,c = 0;

		householder(u,v,b,c,0);

		/* If b = 0, no Householder transformation is necessary */
		if (b != (T) 0)
		{
			mat< T > B(M, M+1, (T) 0);
			mat< T > C(M, m, (T) 0);

			for (size_t i = 0; i < M; i++)
			{
				/* B <-- (I - bvv^t)H(p+1:m-1,p:m-1) */
				for (size_t j = 0; j <= M; j++)
				{
					B(i,j) = H(i+p+1,j+p);

					for (size_t k = 0; k < M; k++)
					{
						B(i,j) -= b * v[i] * v[k] * H(k+p+1,j+p);
					}
				}

				/* C <-- (I - bvv^t)P(p+1:m-1,0:m-1) */
				for (size_t j = 0; j < m; j++)
				{
					C(i,j) = P(i+p+1,j);

					for (size_t k = 0; k < M; k++)
					{
						C(i,j) -= b * v[i] * v[k] * P(k+p+1,j);
					}
				}
			}

			for (size_t i = 0; i < M; i++)
			{
				/* H(p+1:m-1,p:m-1) <--- -c*B */
				for (size_t j = 0; j <= M; j++)
				{
					H(i+p+1,j+p) = -c * B(i,j);
				}

				/* P(p+1:m-1,0:m-1) <--- -c*C */
				for (size_t j = 0; j < m; j++)
				{
					P(i+p+1,j) = -c * C(i,j);
				}
			}

			mat< T > D(m, M, (T) 0);

			/* D <-- H(0:m-1,p+1:m-1)(I - b*vv^t) */
			for (size_t i = 0; i < m; i++)
			{
				for (size_t j = 0; j < M; j++)
				{
					D(i,j) = H(i,j+p+1);

					for (size_t k = 0; k < M; k++)
					{
						D(i,j) -= b * H(i,k+p+1) * v[k] * v[j];
					}
				}
			}

			/* H(0:m-1,p+1:m-1) <--- -c*D */
			for (size_t i = 0; i < m; i++)
			{
				for (size_t j = 0; j < M; j++)
				{
					H(i,j+p+1) = -c*D(i,j);
				}
			}
		}
	}

	/* Q <-- P^t */
	for (size_t i = 0; i < Q.num_rows(); i++)
	{
		for (size_t j = 0; j < Q.num_cols(); j++)
		{
			Q(i,j) = P(j,i);
		}
	}

	return true;
}


template < class T >
template < class _mat1, class _mat2 >
bool hessenberg_house< T >::factorize (const _mat1& A, _mat2& H) const
{
	FREEAML_ASSERT(A.is_square());

	size_t m = A.num_rows();

	/* assign the correct dimensions to H */
	H.resize(m,m);

	/* H <-- A */
	for (size_t i = 0; i < A.num_rows(); i++)
	{
		for (size_t j = 0; j < A.num_cols(); j++)
		{
			H(i,j) = A(i,j);
		}
	}

	/* for each column p of H which must have the portion below
	 * its subdiagonal zeroed out */
	for (size_t p = 0; p < m-1; p++)
	{
		size_t M = m-p-1;

		/* let u := H(p+1:m-1,p) */
		vec< T > u(M, (T) 0);

		for (size_t i = 0; i < M; i++)
		{
			u[i] = H(i+p+1,p);
		}

		/* define the Householder transformation for u which
		 * projects it onto e_0 = (1,0,...0) */
		vec< T > v;

		T b,c = 0;

		householder(u,v,b,c,0);

		/* If b = 0, no Householder transformation is necessary */
		if (b != (T) 0)
		{
			mat< T > B(M, M+1, (T) 0);

			for (size_t i = 0; i < M; i++)
			{
				/* B <-- (I - bvv^t)H(p+1:m-1,p:m-1) */
				for (size_t j = 0; j <= M; j++)
				{
					B(i,j) = H(i+p+1,j+p);

					for (size_t k = 0; k < M; k++)
					{
						B(i,j) -= b * v[i] * v[k] * H(k+p+1,j+p);
					}
				}
			}

			for (size_t i = 0; i < M; i++)
			{
				/* H(p+1:m-1,p:m-1) <--- -c*B */
				for (size_t j = 0; j <= M; j++)
				{
					H(i+p+1,j+p) = -c * B(i,j);
				}
			}

			mat< T > D(m, M, (T) 0);

			/* D <-- H(0:m-1,p+1:m-1)(I - b*vv^t) */
			for (size_t i = 0; i < m; i++)
			{
				for (size_t j = 0; j < M; j++)
				{
					D(i,j) = H(i,j+p+1);

					for (size_t k = 0; k < M; k++)
					{
						D(i,j) -= b * H(i,k+p+1) * v[k] * v[j];
					}
				}
			}

			/* H(0:m-1,p+1:m-1) <--- -c*D */
			for (size_t i = 0; i < m; i++)
			{
				for (size_t j = 0; j < M; j++)
				{
					H(i,j+p+1) = -c*D(i,j);
				}
			}
		}
	}

	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_hessenberg_house_h_ */
