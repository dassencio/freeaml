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


#ifndef _freeAML_rand_h_
#define _freeAML_rand_h_


#include <random>
#include <limits>

#include "debug.h"


#define FREEAML_RANDMIN std::numeric_limits< size_t >::min()
#define FREEAML_RANDMAX std::numeric_limits< size_t >::max()


namespace aml
{

/* pseudorandom number generator */
std::random_device __random_dev;
std:: mt19937_64   __prgen(__random_dev());
std::uniform_int_distribution< size_t > __prdist(FREEAML_RANDMIN, FREEAML_RANDMAX);


/**
 * @brief generates a pseudorandom number on the interval [a,b]
 * @note if a and b are omitted, the random number will be either 0 or 1
 * @note this works for integer types, for other types you should use the
 *       double specialization below
 */
template < class T >
T rand (const T a = (T) 0, const T b = (T) 1)
{
	FREEAML_ASSERT(b > a);

	return a + (T)(__prdist(__prgen) % (size_t)(b - a + 1));
}


/**
 * @brief generates a pseudorandom double on the interval [a,b]
 * @note if a and b are omitted, the random number will be on the range [0,1]
 */
template <>
double rand (const double a, const double b)
{
	FREEAML_ASSERT(b > a);

	return a + (((double) __prdist(__prgen) / (double) FREEAML_RANDMAX) * (b - a));
}


/**
 * @brief generates a pseudorandom float on the interval [a,b]
 * @note if a and b are omitted, the random number will be on the range [0,1]
 */
template <>
float rand (const float a, const float b)
{
	FREEAML_ASSERT(b > a);

	return a + (((float) __prdist(__prgen) / (float) FREEAML_RANDMAX) * (b - a));
}


} /* end of namespace aml */

#endif /* _freeAML_rand_h_ */
