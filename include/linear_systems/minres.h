#ifndef _freeAML_minres_h_
#define _freeAML_minres_h_


#include <general/math.h>

#include "generic_iter.h"


namespace aml
{


template < class T >
class minres: public generic_iter< T >
{

public:

	typedef generic_iter< T >	gen_it_type;


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
	minres (const size_t max_iter, const T& residual_tol);


	/***********************************************************************
	 *
	 *	SOLVE FUNCTION
	 *
	 **********************************************************************/

	/**
	 * @brief solves a linear system using the MINRES method
	 * @param A the system matrix
	 * @param x a vector on which the solution will be written (if found)
	 * @param b the right-hand side vector
	 * @return true if the residual tolerance could be achieved within the
	 *         maximum number of iterations allowed, false otherwise
	 */
	template < class _mat, class _vec >
	bool solve (const _mat& A, _vec& x, const _vec& b);

}; /* end of class minres */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
minres< T >::minres(const size_t max_iter, const T& residual_tol):

        gen_it_type(max_iter, residual_tol)
{
	/* nothing needs to  be done here */
}


template < class T >
template < class _mat, class _vec >
bool minres< T >::solve (const _mat& A, _vec& x, const _vec& b)
{
	FREEAML_ASSERT(gen_it_type::check_dimensions(A,x,b));

	size_t n = A.num_rows();

	/* in this algorithm, the initial guess is the zero vector */
	x.zero_fill();

	/*  q_{k-1}, q_k and q_{k+1}, k is the iteration number */
	_vec q1(n);
	_vec q2(n);
	_vec q3(n);

	q1.zero_fill();
	q2.zero_fill();
	q3.zero_fill();

	/* m_{k-2}, m_{k-1} and m_k */
	_vec m1(n);
	_vec m2(n);
	_vec m3(n);

	/* alpha_k, beta_k, beta_{k+1} (these represent the only nonzero
	 * elements of the last column of H_k) */
	T a1 = (T) 0;
	T b1 = (T) 0;
	T b2 = (T) 0;

	/* technically we have here k = -1, q_{k+1} = q_0 = b is the first
	 * column of the Kryvlov subspace basis (matrix) Q */
	q3 = b;

	b1 = q3.l2_norm();

	/* if the right-hand side b is zero, the solution is then x = 0 */
	if (b1 == (T) 0)
	{
		gen_it_type::set_num_iterations(0);
		return true;
	}

	/* normalize q3 */
	q3 /= b1;

	/* t_k, the k-th component of the vector t of the MINRES algorithm */
	T t1 = (T) b1;

	/* Givens rotation parameters which define G(k-2,k-1) */
	T c1 = (T) 0;
	T s1 = (T) 0;

	/* Givens rotation parameters which define G(k-1,k) */
	T c2 = (T) 0;
	T s2 = (T) 0;

	/* Givens rotation parameters which define G(k,k+1) */
	T c3 = (T) 0;
	T s3 = (T) 0;

	size_t num_iter = 0;

	while (num_iter <= gen_it_type::max_num_iterations())
	{
		size_t k = num_iter;

		/* q_{k-1} <--- q_k and q_k <--- q_{k+1} */
		q1 = q2;
		q2 = q3;

		/* q_{k+1} <--- A*q_k */
		q3 = A*q2;

		/* alpha_k <--- (q_k)^t*q_{k+1} */
		a1 = q2*q3;

		/* q_{k+1} <--- q_{k+1} - alpha_k*q_k - beta_k*q_{k-1} */
		q3 -= a1*q2 + b1*q1;

		/* beta_{k+1} <--- norm(q_{k+1})_2; last column of H_k is
		 * beta_k, alpha_k, beta_{k+1} */
		b2 = q3.l2_norm();

		q3 /= b2;

		/* let (epsilon_k,delta_k,gamma_k) = (0,beta_k,alpha_k) */
		T ek = (T) 0;
		T dk = b1;
		T gk = a1;

		/* make copy zeta_k of beta_{k+1} */
		T zk = b2;

		m3 = q2;

		if (k > 1)
		{
			/* rotate (epsilon_k,delta_k) using G(k-2,k-1) to obtain
			 * epsilon_k; delta_k is temporary */
			givens_rotation(ek, dk, c1, s1);

			m3 -= ek*m1;
		}
		if (k > 0)
		{
			/* rotate (delta_k,gamma_k) using G(k-1,k) to obtain
			 * delta_k; gamma_k is temporary */
			givens_rotation(dk, gk, c2, s2);

			m3 -= dk*m2;
		}

		/* get the Givens rotation matrix parameters to make the term
		 * T(k+1,k) equal to zero */
		givens(gk, b2, c3, s3);

		/* rotate (gamma_k,zeta_k) using G(k,k+1) to obtain gamma_k
		 * and zero out beta_{k+1} on T */
		givens_rotation(gk, zk, c3, s3);

		m3 /= gk;

		/* t_{k+1}, the (k+1)-th component of the vector t of the
		 * MINRES algorithm */
		T t2 = (T) 0;

		/*
		 * apply the Givens rotation G(k,k+1) to the last two nonzero
		 * elements of the right-hand side: t_{k+1} <--- t_k*s_k
		 * and t_k <--- t_k*c_k
		 */
		givens_rotation(t1, t2, c3, s3);

		x += t1*m3;

		/* t_k <--- t_{k+1} */
		t1 = t2;

		/* m_{k-2} <--- m_{k-1} and m_{k-1} <--- m_k */
		m1 = m2;
		m2 = m3;

		c1 = c2;
		s1 = s2;

		c2 = c3;
		s2 = s3;

		/* beta_k <--- beta_{k+1} */
		b1 = b2;

		num_iter++;

		/* if the residual is within the maximum tolerance, stop */
		if (std::abs(t1) <= gen_it_type::residual_tolerance())
		{
			gen_it_type::set_num_iterations(num_iter);
			return true;
		}
	}

	gen_it_type::set_num_iterations(num_iter);
	return false;
}


} /* end of namespace aml */

#endif /* _freeAML_minres_h_ */
