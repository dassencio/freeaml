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


#ifndef _freeAML_math_h_
#define _freeAML_math_h_


#include <cstdlib>

#include "debug.h"


namespace aml
{


/*******************************************************************************
 *
 *	BASIC MATHEMATICAL FUNCTIONS
 *
 ******************************************************************************/

/**
 * @brief computes the Kronecker delta for a pair of integers
 * @param i an integer
 * @param j an integer
 * @return one if i and j are equal or zero otherwise
 */
int kronecker_delta (const int i, const int j)
{
	return (i == j) ? 1 : 0;
}


/**
 * @brief checks whether an integer is a power of 2
 * @param n an integer
 * @return true if n is a power of 2, false otherwise
 */
inline bool is_power_of_2 (const int n)
{
	return (n > 0) && !(n & (n - 1));
}


/**
 * @brief checks whether an integer is even
 * @param n an integer
 * @return true if n is even, false otherwise
 */
inline bool is_even (const int n)
{
	return (n % 2 == 0);
}


/**
 * @brief checks whether an integer is odd
 * @param n an integer
 * @return true if n is odd, false otherwise
 */
inline bool is_odd (const int n)
{
	return !is_even(n);
}


/**
 * @brief computes the sign of a number
 * @param x a number (can be int, double, float etc.)
 * @return 1 if x is nonnegative, -1 otherwise
 */
template < class T >
inline T sign (const T& x)
{
	if (x >= 0) return (T) 1;
	return (T) -1;
}


/**
 * @brief computes the number of digits of an input integer
 * @param n an integer
 * @return the number of digits on n (if n is negative, the minus sign is ignored)
 */
inline size_t num_digits (int n)
{
	/* every integer has at least one digit */
	int digits = 1;

	/* as long as n has more than one digit */
	while (n / 10 != 0)
	{
		digits ++;

		/* throw away its last digit */
		n /= 10;
	}
	return digits;
}


/**
 * @brief computes a binomial coefficient
 * @param n an integer
 * @param k an integer
 * @return the binomial coefficient indexed by n and k: n!/(k!(n-k)!) if
 *         0 <= k <= n, otherwise zero
 */
int binomial (const int n, const int k)
{
	/* hard-coded cases: n = 0,...6; 0 <= k <= n */
	switch (n)
	{
	case 0:
		switch (k)
		{
		case 0:  return 1;
		default: return 0;
		}
	case 1:
		switch (k)
		{
		case 0:  return 1;
		case 1:  return 1;
		default: return 0;
		}
	case 2:
		switch (k)
		{
		case 0:  return 1;
		case 1:  return 2;
		case 2:  return 1;
		default: return 0;
		}
	case 3:
		switch (k)
		{
		case 0:  return 1;
		case 1:  return 3;
		case 2:  return 3;
		case 3:  return 1;
		default: return 0;
		}
	case 4:
		switch (k)
		{
		case 0:  return 1;
		case 1:  return 4;
		case 2:  return 6;
		case 3:  return 4;
		case 4:  return 1;
		default: return 0;
		}
	case 5:
		switch (k)
		{
		case 0:  return 1;
		case 1:  return 5;
		case 2:  return 10;
		case 3:  return 10;
		case 4:  return 5;
		case 5:  return 1;
		default: return 0;
		}
	case 6:
		switch (k)
		{
		case 0:  return 1;
		case 1:  return 6;
		case 2:  return 15;
		case 3:  return 20;
		case 4:  return 15;
		case 5:  return 6;
		case 6:  return 1;
		default: return 0;
		}
	default:
		break;
	}

	if(k < 0 || k > n)
	{
		return 0;
	}
	else
	{
		/* compute all other cases by hand */
		int rv = 1;

		for (int i = 1; i <= k; i++)
		{
			rv *= (n - k + i) / i;
		}

		return rv;
	}
}


/*******************************************************************************
 *
 *	OTHER MATHEMATICAL FUNCTIONS
 *
 ******************************************************************************/

/**
 * @brief computes the Householder vector and the associated parameters which
 *        define the Householder transformation for a given input vector
 * @param u the vector for which the Householder transformation will be defined
 * @param v at the end, the Householder transformation vector for u
 * @param b at the end, the first Householder transformation constant
 * @param c at the end, the second Householder transformation constant
 * @param j the index such that the Householder transformation for u projects
 *        it onto the unit vector e_j whose elements are all zero except for
 *        the j-th element which is 1 (if omitted, j is assumed to be zero)
 * @note the vector v and the constants b and c are such that
 *       (I - bvv^t)u = -c|u|e_j, where I is the identity matrix, v is a column
 *       vector and |u| is the L^2 norm of the vector u
 */
template < class T, class _vec >
void householder (const _vec& u, _vec& v, T& b, T& c, const size_t j = 0)
{
	v = u;

	FREEAML_ASSERT(j < v.size());

	T norm = u.l2_norm();

	/* if u is the null vector, then it is already a multiple of e_j */
	if (norm == (T) 0)
	{
		b = (T) 0;
		c = (T) 1;
	}
	else
	{
		c = sign(u[j]);

		/* only the j-th element of v needs to be adjusted */
		v[j] += c * norm;

		v /= v.l2_norm();

		b = (T) 2;
	}
}


/**
 * @brief computes the entries of the Givens rotation matrix for a 2d vector
 * @param x the x component of the 2d vector
 * @param y the y component of the 2d vector
 * @param c the first parameter which defines the Given rotation matrix for (x,y)
 * @param s the second parameter which defines the Given rotation matrix for (x,y)
 * @note the Givens rotation matrix has the form G = ( c -s ; s c )
 */
template < class T >
void givens (const T& x, const T& y, T& c, T& s)
{
	/* compute the magnitude of (x,y) */
	T norm = std::sqrt(x*x + y*y);

	if (norm > (T) 0)
	{
		c =  x / norm;
		s = -y / norm;
	}
	else
	{
		c = (T) 0;
		s = (T) 0;
	}
}


} /* end of namespace aml */

#endif /* _freeAML_math_h_ */
