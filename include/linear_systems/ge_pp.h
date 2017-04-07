#ifndef _freeAML_ge_pp_h_
#define _freeAML_ge_pp_h_


#include <cmath>
#include <algorithm>

#include <general/debug.h>


/*******************************************************************************
 *
 *	GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING (GEPP)
 *
 ******************************************************************************/


namespace aml
{


template < class T >
class ge_pp
{

public:

	/**
	 * @brief solves a linear system using gaussian elimination
	 *        with partial pivoting
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @return true if the linear system could be solved, false otherwise
	 */
	template < class _mat, class _vec >
	bool solve (_mat A, _vec& x, _vec b) const;

}; /* end of class ge_pp */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
template < class _mat, class _vec >
bool ge_pp< T >::solve (_mat A, _vec& x, _vec b) const
{
	FREEAML_ASSERT(A.is_square());
	FREEAML_ASSERT(x.size() == b.size());
	FREEAML_ASSERT(A.num_rows() == x.size());

	size_t n = A.num_rows();

	_vec row(n);

	/* row index vector: used for "swapping" rows of A */
	for (size_t i = 0; i < n; i++)
	{
		row[i] = i;
	}

	/* for each row i of the matrix A */
	for (size_t i = 0; i < n; i++)
	{
		size_t p = i;
		size_t q = i;

		/* find p such that |A_pi| = max_q|A_row[q],i)|, i <= q < n */
		while (q < n)
		{
			if (std::abs(A(row[p],i)) < std::abs(A(row[q],i)))
			{
				p = q;
			}
			q++;
		}

		/* if A_(row[q],i) = 0, A is not invertible */
		if (A(row[p],i) == (T) 0)
		{
			return false;
		}

		if (p != i)
		{
			/* "swap" rows i (A_i) and j (A_j) of A */
			std::swap(row[p],row[i]);
		}

		for (size_t j = i+1; j < n; j++)
		{
			T m = A(row[j],i) / A(row[i],i);

			/* A_j <--- A_j - m*A_i */
			for (size_t l = i+1; l < n; l++)
			{
				A(row[j],l) -= m * A(row[i],l);
			}

			b[row[j]] -= m*b[row[i]];
		}
	}

	/* backward substitution */
	size_t i = n;

	while (i-- > 0)
	{
		T y = (T) 0;

		for (size_t j = i+1; j < n; j++)
		{
			y += A(row[i],j) * x[j];
		}
		x[i] = (b[row[i]] - y) / A(row[i],i);
	}

	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_ge_pp_h_ */
