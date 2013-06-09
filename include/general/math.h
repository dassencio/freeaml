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


#ifndef _freeAML_basic_math_h_
#define _freeAML_basic_math_h_


namespace aml
{


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
	FREEAML_ASSERT(k >= 0);
	FREEAML_ASSERT(n >= k);

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

	/* compute all other cases by hand */
	int rv = 1;

	for (int i = 1; i <= k; i++)
	{
		rv *= (n - k + i) / i;
	}
}




} /* end of namespace aml */

#endif /* _freeAML_basic_math_h_ */
