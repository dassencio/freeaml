#ifndef _freeAML_cholesky_h_
#define _freeAML_cholesky_h_


#include <cmath>

#include <general/debug.h>


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
