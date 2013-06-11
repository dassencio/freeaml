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


#ifndef _freeAML_qr_h_
#define _freeAML_qr_h_


#include <cstdlib>

#include "data_structures/mat.h"


namespace aml
{


template < class T >
class qr
{

public:

	/**
	 * @brief solves a least squares problem using a QR factorization method
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @param fac a QR factorization class
	 * @return true if a solution could be found, false otherwise
	 */
	template < class _mat, class _vec, class _qr >
	bool solve (const _mat& A, _vec& x, const _vec& b, const _qr& fac);

}; /* end of class qr */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
template < class _mat, class _vec, class _qr >
bool qr< T >::solve (const _mat& A, _vec& x, const _vec&  b, const _qr& fac)
{
	/* A is an m × n matrix */
	size_t n = A.num_cols();

	/* A must be either nonsingular or represent an overdetermined
	 * linear system */
	FREEAML_ASSERT(A.num_rows() >= n);

	/* define the matrices Q and R (they will be resized automatically
	 * by the QR factorization class) */
	aml::mat< T > Q,R;

	/* compute the QR factorization A = QR */
	if (!fac.factorize(A,Q,R))
	{
		return false;
	}

	/* let f = (Q^t)*b */
	_vec f = Q.transpose() * b;

	/* solve R*x = b using backward substitution (R should be upper
	 * triangular with nonzero diagonal elements) */
	size_t i = n;

	while (i-- > 0)
	{
		T y = (T) 0;

		for (size_t j = i+1; j < n; j++)
		{
			y += R(i,j) * x[j];
		}

		if (R(i,i) == (T) 0)
		{
			return false;
		}

		x[i] = (f[i] - y) / R(i,i);
	}

	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_qr_h_ */
