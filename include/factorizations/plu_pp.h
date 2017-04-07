#ifndef _freeAML_plu_pp_h_
#define _freeAML_plu_pp_h_


#include <data_structures/vec.h>


namespace aml
{


template < class T >
class plu_pp
{

public:

	/***********************************************************************
	 *
	 *	FACTORIZE FUNCTIONS
	 *
	 **********************************************************************/

	/**
	 * @brief computes the PLU factorization of an input matrix with
	 *        partial pivoting
	 * @param A the matrix which must be factorized
	 * @param P on success, a permutation matrix
	 * @param L on success, a lower triangular matrix with all diagonal
	 *        elements set to one (1)
	 * @param U on success, an upper triangular nonsingular matrix
	 * @return true if the factorization A = PLU could be built, false otherwise
	 * @note the PLU factorization exists only if A is invertible
	 * @note partial pivoting is used when building the PLU factorization
	 */
	template < class _mat1, class _mat2,
	           class _mat3, class _mat4 >
	bool factorize (_mat1  A,
	                _mat2& P,
	                _mat3& L,
	                _mat4& U) const;

}; /* end of class plu */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
template < class _mat1, class _mat2,
           class _mat3, class _mat4 >
bool plu_pp< T >::factorize (_mat1  A,
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

	/* perm will be used to build the permutation matrix P at the end */
	vec< size_t > perm(n,0);

	/* initialize perm as the "indentity permutation" */
	for (size_t p = 0; p < n; p++)
	{
		perm[p] = p;
	}

	/* for each row p of A (excluding for the last one) */
	for (size_t p = 0; p < n-1; p++)
	{
		/* do partial pivoting to make A(perm[p],p) as large (in
		 * magnitude) as possible */

		for (size_t j = p+1; j < n; j++)
		{
			if (std::abs(A(perm[j],p)) > std::abs(A(perm[p],p)))
			{
				std::swap(perm[p],perm[j]);
			}
		}

		/* if no nonzero entry was found on or below the diagonal
		 * element of the p-th row of A, then A is not invertible */
		if (A(perm[p],p) == (T) 0)
		{
			return false;
		}

		for (size_t k = p+1; k < n; k++)
		{
			A(perm[k],p) /= A(perm[p],p);

			for (size_t i = p+1; i < n; i++)
			{
				A(perm[k],i) -= A(perm[k],p) * A(perm[p],i);
			}
		}
	}

	/* define the permutation matrix P and the diagonal of L */
	for (size_t i = 0; i < n; i++)
	{
		L(i,i) = (T) 1;

		P(perm[i],i) = (T) 1;
	}

	/* define the rest of L */
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < i; j++)
		{
			L(i,j) = A(perm[i],j);
		}
	}

	/* define U */
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = i; j < n; j++)
		{
			U(i,j) = A(perm[i],j);
		}
	}

	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_plu_pp_h_ */
