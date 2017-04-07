#ifndef _freeAML_plu_h_
#define _freeAML_plu_h_


#include <data_structures/vec.h>


namespace aml
{


template < class T >
class plu
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
bool plu< T >::factorize (_mat1  A,
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
		/* find the first nonzero entry on or below the diagonal
		 * element of the p-th row of A */

		size_t i = p;

		while (i < n && A(perm[i],p) == (T) 0)
		{
			i++;
		}

		/* if no such i exists, A is not invertible */
		if (i == n)
		{
			return false;
		}

		if (i != p)
		{
			std::swap(perm[i],perm[p]);
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

#endif /* _freeAML_plu_h_ */
