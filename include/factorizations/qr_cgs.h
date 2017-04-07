#ifndef _freeAML_qr_cgs_h_
#define _freeAML_qr_cgs_h_


#include <cmath>

#include <general/debug.h>


namespace aml
{


template < class T >
class qr_cgs
{

public:

	/**
	 * @brief computes the QR factorization of a matrix using the classical
	 *        version of the Gram-Schmidt method
	 * @param A the matrix which must be factorized
	 * @param Q on success, the matrix Q of the QR factorization
	 * @param R on success, the matrix R of the QR factorization
	 * @return true if the factorization A = QR could be built, false otherwise
	 * @note if the method succeeds, Q will be orthogonal and R will be
	 *       upper triangular with positive diagonal values
	 * @note the number of rows of A must be larger or equal than its number
	 *       of columns
	 */
	template < class _mat1, class _mat2, class _mat3 >
	bool factorize (const _mat1& A, _mat2& Q, _mat3& R) const;

}; /* end of class qr_cgs */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
template < class _mat1, class _mat2, class _mat3 >
bool qr_cgs< T >::factorize (const _mat1& A, _mat2& Q, _mat3& R) const
{
	size_t m = A.num_rows();
	size_t n = A.num_cols();

	/* A must be either nonsingular or represent an overdetermined
	 * linear system */
	FREEAML_ASSERT(m >= n);

	/* assign the correct dimensions to Q and R */
	Q.resize(m,n);
	R.resize(n,n);

	Q.zero_fill();
	R.zero_fill();

	for (size_t i = 0; i < n; i++)
	{
		/* make the i-th column of Q, (q_i) equal to the i-th
		 * column of A (a_i) */
		for (size_t j = 0; j < m; j++)
		{
			Q(j,i) = A(j,i);
		}

		for (size_t j = 0; j < i; j++)
		{
			/* R_{ji} <--- q_j^t * a_i (dot product) */
			for (size_t k = 0; k < m; k++)
			{
				R(j,i) += Q(k,j) * A(k,i);
			}

			/* q_i <--- q_i - R_{ji}*q_j */
			for (size_t k = 0; k < m; k++)
			{
				Q(k,i) -= R(j,i) * Q(k,j);
			}
		}

		/* compute the L^2 norm of q_i */
		T norm = 0;

		for (size_t j = 0; j < m; j++)
		{
			norm += Q(j,i) * Q(j,i);
		}

		norm = std::sqrt(norm);

		R(i,i) = norm;

		if (norm == (T) 0)
		{
			return false;
		}

		/* normalize q_i */
		for (size_t j = 0; j < m; j++)
		{
			Q(j,i) /= norm;
		}
	}

	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_qr_cgs_h_ */
