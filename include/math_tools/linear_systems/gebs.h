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


#ifndef _freeAML_gebs_h_
#define _freeAML_gebs_h_


#include <cmath>
#include <algorithm>

#include "general/debug.h"


/*******************************************************************************
 *
 *	GAUSSIAN ELIMINATION WITH BACKWARD SUBSTITUTION (GEBS)
 *
 ******************************************************************************/


namespace aml
{


template < class T >
class gebs
{

public:

	/**
	 * @brief solves a linear system using regular gaussian
	 *        elimination (no sophisticated pivoting strategy)
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @return true if the linear system could be solved, false otherwise
	 */
	template < class _mat, class _vec >
	bool solve (_mat A, _vec& x, _vec b) const;

}; /* end of class gebs */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
template < class _mat, class _vec >
bool gebs< T >::solve (_mat A, _vec& x, _vec b) const
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

		/* find the first p such that A_row[p],i != 0 */
		while (p < n && A(row[p],i) == (T) 0)
		{
			p++;
		}

		/* if no such p exists, A is not invertible */
		if (p == n)
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

#endif /* _freeAML_gebs_h_ */
