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


#ifndef _freeAML_jp_minres_h_
#define _freeAML_jp_minres_h_


#include "minres.h"

#include <data_structures/sparse_mat.h>


namespace aml
{


template < class T >
class jp_minres: public minres< T >
{

public:

	/***********************************************************************
	 *
	 *	CONSTRUCTORS
	 *
	 **********************************************************************/

	/**
	 * @brief default constructor
	 * @param max_iter maximum number of iterations allowed in one solve
	 * @param residual_tol maximum residual tolerance allowed
	 */
	jp_minres (const size_t max_iter, const T& residual_tol);


	/***********************************************************************
	 *
	 *	SOLVE FUNCTION
	 *
	 **********************************************************************/

	/**
	 * @brief solves a linear system using the Jacobi preconditioned
	 *        MINRES method
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @return true if the residual tolerance could be achieved within the
	 *         maximum number of iterations allowed, false otherwise
	 * @note a Jacobi preconditioner can only be used if all diagonal
	 *       elements of A are nonzero
	 */
	template < class _mat, class _vec >
	bool solve (const _mat& A, _vec& x, const _vec& b);

}; /* end of class jp_minres */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
jp_minres< T >::jp_minres(const size_t max_iter, const T& residual_tol):

        minres< T >(max_iter, residual_tol)
{
	/* nothing needs to  be done here */
}


template < class T >
template < class _mat, class _vec >
bool jp_minres< T >::solve (const _mat& A, _vec& x, const _vec& b)
{
	/* D is the Jacobi preconditioner matrix */
	sparse_mat< T > D(A.num_rows(), A.num_cols());

	minres< T >::set_num_iterations(0);

	/* build D */
	for (size_t i = 0; i < D.num_rows(); i++)
	{
		/* if A(i,i) = 0, the Jacobi preconditioner cannot be built */
		if (A(i,i) == (T) 0)
		{
			return false;
		}
		else
		{
			D(i,i) = (T) 1 / std::sqrt(std::abs(A(i,i)));
		}
	}

	/* solve (DAD)*x = D*b (preconditioned linear system); here x is
	 * not the final solution; D*x is! */
	if (!minres< T >::solve(D*A*D, x, D*b))
	{
		return false;
	}

	x = D*x;

	return true;
}


} /* end of namespace aml */

#endif /* _freeAML_jp_minres_h_ */
