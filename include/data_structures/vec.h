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


#ifndef _freeAML_vec_h_
#define _freeAML_vec_h_


#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>

#include "general/debug.h"
#include "general/math.h"


namespace aml
{


template < class T >
class vec: protected std::vector< T >
{

public:

	typedef std::vector< T >	vector_type;


	/***********************************************************************
	 *
	 *	USED MEMBERS OF BASE CLASS std::vector
	 *
	 **********************************************************************/

	using vector_type::assign;
	using vector_type::at;
	using vector_type::back;
	using vector_type::begin;
	using vector_type::capacity;
	using vector_type::clear;
	using vector_type::const_iterator;
	using vector_type::const_pointer;
	using vector_type::const_reference;
	using vector_type::const_reverse_iterator;
	using vector_type::data;
	using vector_type::difference_type;
	using vector_type::empty;
	using vector_type::end;
	using vector_type::erase;
	using vector_type::front;
	using vector_type::get_allocator;
	using vector_type::insert;
	using vector_type::iterator;
	using vector_type::max_size;
	using vector_type::pointer;
	using vector_type::pop_back;
	using vector_type::push_back;
	using vector_type::rbegin;
	using vector_type::reference;
	using vector_type::rend;
	using vector_type::reserve;
	using vector_type::resize;
	using vector_type::reverse_iterator;
	using vector_type::size;
	using vector_type::size_type;
	using vector_type::value_type;


	/***********************************************************************
	 *
	 *	CONSTRUCTORS
	 *
	 **********************************************************************/

	/** @brief default constructor */
	explicit vec ();

	/** @brief fill constructor */
	explicit vec (size_t len, const T& defval = T());

	/** @brief range constructor */
	template < class input_iterator >
	vec (input_iterator first, input_iterator last);

	/** @brief copy constructor */
	vec (const vec< T >& v);


	/***********************************************************************
	 *
	 *	OVERLOADED OPERATORS
	 *
	 **********************************************************************/

	/** @brief gets a reference to the i-th element of the vector */
	T& operator[] (const size_t i);

	/** @brief gets a const reference to the i-th element of the vector */
	const T& operator[] (const size_t i) const;

	/** @brief multiplies all elements of the vector by a scalar */
	vec< T >& operator*= (const T& c);

	/** @brief computes the multiplication of the vector by a scalar (on the right) */
	vec< T > operator* (const T& c) const;

	/** @brief divides all elements of the vector by a scalar */
	vec< T >& operator/= (const T& c);

	/** @brief computes the division of the vector by a scalar */
	vec< T > operator/ (const T& c) const;

	/** @brief adds an input vector to this vector */
	vec< T >& operator+= (const vec< T >& v);

	/** @brief computes the addition of this vector and an input vector */
	vec< T > operator+ (const vec< T >& v) const;

	/** @brief subtracts an input vector from this vector */
	vec< T >& operator-= (const vec< T >& v);

	/** @brief computes the subtraction of an input vector from this vector */
	vec< T > operator- (const vec< T >& v) const;

	/** @brief computes the dot product of this vector and an input vector */
	T operator* (const vec< T >& v) const;


	/***********************************************************************
	 *
	 *	MATHEMATICAL FUNCTIONS
	 *
	 **********************************************************************/

	/** @brief computes the L^1 norm of the vector */
	T l1_norm () const;

	/** @brief computes the L^2 norm of the vector */
	T l2_norm () const;

	/** @brief computes the L^p norm of the vector */
	T lp_norm (const int p) const;

	/** @brief computes the L^inf norm of the vector */
	T linf_norm () const;

	/** @brief computes the sum of all elements of the vector */
	T sum (void) const;

	/** @brief computes the average of all elements of the vector */
	T average (void) const;

	/** @brief changes the vector average value to an input value */
	vec< T >& average_elements (const T& x);


	/***********************************************************************
	 *
	 *	OTHER FUNCTIONS
	 *
	 **********************************************************************/

	/** @brief sets all elements of the vector to a constant */
	void fill (const T& x);

	/** @brief sets all elements of the vector to zero */
	void zero_fill (void);

	/** @brief swaps the contents of this vector and an input vector */
	void swap (vec< T >& v);


	/***********************************************************************
	 *
	 *	PRINT FUNCTIONS
	 *
	 **********************************************************************/

	/**
	 * @brief prints all elements of the vector
	 * @param st the (optional) output stream
	 * @return the output stream st
	 */
	std::ostream& print (std::ostream& st = std::cout) const;

	/**
	 * @brief prints all elements of the vector with their indices
	 * @param st the (optional) output stream
	 * @return the output stream st
	 */
	std::ostream& print_with_indices (std::ostream& st = std::cout) const;

}; /* end of class vec */


/*******************************************************************************
 *
 *	NON-MEMBER OVERLOADED OPERATORS
 *
 ******************************************************************************/

/** @brief computes the multiplication of a vector by a scalar */
template < class T >
vec< T > operator* (const T& c, const vec< T >& v);

/** @brief prints the elements of a vector directly to an output stream */
template < class T >
std::ostream& operator<< (std::ostream& st, const vec< T >& v);


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS AND MEMBER OPERATORS
 *
 ******************************************************************************/

template < class T >
vec< T >::vec (): vector_type()
{
	/* nothing needs to be done here */
}


template < class T >
vec< T >::vec (size_t len, const T& defval): vector_type(len,defval)
{
	/* nothing needs to be done here */
}


template < class T >
template < class input_iterator >
vec< T >::vec (input_iterator first, input_iterator last): vector_type(first,last)
{
	/* nothing needs to be done here */
}


template < class T >
vec< T >::vec (const vec< T > &v): vector_type(v.begin(),v.end())
{
	/* nothing needs to be done here */
}


template < class T >
T& vec< T >::operator[] (const size_t i)
{
	FREEAML_ASSERT(i < size());

	return vector_type::operator[](i);
}


template < class T >
const T& vec< T >::operator[] (const size_t i) const
{
	FREEAML_ASSERT(i < size());

	return vector_type::operator[](i);
}


template < class T >
vec< T >& vec< T >::operator*= (const T& c)
{

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < size(); i++)
	{
		(*this)[i] *= c;
	}
	return (*this);
}


template < class T >
vec< T > vec< T >::operator* (const T& c) const
{
	vec< T > u = (*this);
	return (u *= c);
}


template < class T >
vec< T >& vec< T >::operator/= (const T& c)
{

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < size(); i++)
	{
		(*this)[i] /= c;
	}
	return (*this);
}


template < class T >
vec< T > vec< T >::operator/ (const T& c) const
{
	vec< T > u = (*this);
	return (u /= c);
}


template < class T >
vec< T >& vec< T >::operator+= (const vec< T >& v)
{
	FREEAML_ASSERT(size() == v.size());

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < size(); i++)
	{
		(*this)[i] += v[i];
	}
	return (*this);
}


template < class T >
vec< T > vec< T >::operator+ (const vec< T >& v) const
{
	FREEAML_ASSERT(size() == v.size());

	vec< T > u = (*this);
	return (u += v);
}


template < class T >
vec< T >& vec< T >::operator-= (const vec< T >& v)
{
	FREEAML_ASSERT(size() == v.size());

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < size(); i++)
	{
		(*this)[i] -= v[i];
	}
	return (*this);
}


template < class T >
vec< T > vec< T >::operator- (const vec< T >& v) const
{
	FREEAML_ASSERT(size() == v.size());

	vec< T > u = (*this);
	return (u -= v);
}


template < class T >
T vec< T >::operator* (const vec< T >& v) const
{
	FREEAML_ASSERT(size() == v.size());

	T dotp = (T) 0;

#ifdef _OPENMP
#pragma omp parallel
	{
		T local_dotp = (T) 0;

#pragma omp for nowait

		for (size_t i = 0; i < size(); i++)
		{
			local_dotp += (*this)[i] * v[i];
		}

#pragma omp critical
		{
			dotp += local_dotp;
		}
	}
#else

	/* serial implementation */
	for (size_t i = 0; i < size(); i++)
	{
		dotp += (*this)[i] * v[i];
	}
#endif

	return dotp;
}


template < class T >
T vec< T >::l1_norm () const
{
	T norm = (T) 0;

#ifdef _OPENMP
#pragma omp parallel
	{
		T local_norm = (T) 0;

#pragma omp for nowait

		for (size_t i = 0; i < size(); i++)
		{
			local_norm += std::abs((*this)[i]);
		}

#pragma omp critical
		{
			norm += local_norm;
		}
	}
#else

	/* serial implementation */
	for (size_t i = 0; i < size(); i++)
	{
		norm += std::abs((*this)[i]);
	}
#endif

	return norm;
}


template < class T >
T vec< T >::l2_norm () const
{
	T norm = (T) 0;

#ifdef _OPENMP
#pragma omp parallel
	{
		T local_norm = (T) 0;

#pragma omp for nowait

		for (size_t i = 0; i < size(); i++)
		{
			local_norm += std::abs((*this)[i]) * std::abs((*this)[i]);
		}

#pragma omp critical
		{
			norm += local_norm;
		}
	}
#else

	/* serial implementation */
	for (size_t i = 0; i < size(); i++)
	{
		norm += std::abs((*this)[i]) * std::abs((*this)[i]);
	}
#endif

	return std::sqrt(norm);
}


template < class T >
T vec< T >::lp_norm (const int p) const
{
	T norm = (T) 0;

#ifdef _OPENMP
#pragma omp parallel
	{
		T local_norm = (T) 0;

#pragma omp for nowait

		for (size_t i = 0; i < size(); i++)
		{
			local_norm += std::pow(std::abs((*this)[i]), p);
		}

#pragma omp critical
		{
			norm += local_norm;
		}
	}
#else

	/* serial implementation */
	for (size_t i = 0; i < size(); i++)
	{
		norm += std::pow(std::abs((*this)[i]), p);
	}
#endif

	return std::pow(norm, (T)1/p);
}


template < class T >
T vec< T >::linf_norm () const
{
	T norm = (T) 0;

#ifdef _OPENMP
#pragma omp parallel
	{
		T local_norm = (T) 0;

#pragma omp for nowait

		for (size_t i = 0; i < size(); i++)
		{
			/* this should also work with complex numbers */
			local_norm = std::max(std::abs(local_norm), std::abs((*this)[i]));
		}

#pragma omp critical
		{
			norm = std::max(std::abs(local_norm), std::abs(norm));
		}
	}
#else

	/* serial implementation */
	for (size_t i = 0; i < size(); i++)
	{
		norm = std::max(std::abs(norm), std::abs((*this)[i]));
	}
#endif

	return norm;
}


template < class T >
T vec< T >::sum (void) const
{
	T sum = (T) 0;

#ifdef _OPENMP
#pragma omp parallel
	{
		T local_sum = (T) 0;

#pragma omp for nowait

		for (size_t i = 0; i < size(); i++)
		{
			local_sum += (*this)[i];
		}

#pragma omp critical
		{
			sum += local_sum;
		}
	}
#else

	/* serial implementation */
	for (size_t i = 0; i < size(); i++)
	{
		sum += (*this)[i];
	}
#endif

	return sum;
}


template < class T >
T vec< T >::average (void) const
{
	if (!empty())
	{
		return sum() / (T) size();
	}
	else
	{
		return (T) 0;
	}
}


template < class T >
vec< T >& vec< T >::average_elements (const T& x)
{
	const T avg = average();

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < size(); i++)
	{
		(*this)[i] += (x - avg);
	}
	return *this;
}


template < class T >
void vec< T >::fill (const T& x)
{

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < size(); i++)
	{
		(*this)[i] = x;
	}
}


template < class T >
void vec< T >::zero_fill (void)
{
	fill((T)0);
}


template < class T >
void vec< T >::swap (vec< T >& v)
{
	vector_type::swap(v);
}


template < class T >
std::ostream& vec< T >::print (std::ostream& st) const
{
	size_t max_width = 0;

	for (size_t i = 0; i < size(); i++)
	{
		std::stringstream ss;
		ss << (*this)[i];

		max_width = std::max(max_width, ss.str().size());
	}

	for (size_t i = 0; i < size(); i++)
	{
		st << std::setw(max_width+2) << (*this)[i] << "\n";
	}

	return st;
}


template < class T >
std::ostream& vec< T >::print_with_indices (std::ostream& st) const
{
	size_t max_width = 0;

	/* get the number of digits of the vector size */
	size_t n = num_digits(size());

	for (size_t i = 0; i < size(); i++)
	{
		std::stringstream ss;
		ss << (*this)[i];

		max_width = std::max(max_width, ss.str().size());
	}

	for (size_t i = 0; i < size(); i++)
	{
		st << "["
		   << std::setw(n) << i
		   << "] "
		   << std::setw(max_width+2)
		   << (*this)[i] << "\n";
	}

	return st;
}


/*******************************************************************************
 *
 *	DEFINITIONS OF NON-MEMBER OVERLOADED OPERATORS
 *
 ******************************************************************************/

template < class T >
vec< T > operator* (const T& c, const vec< T >& v)
{
	return (v*c);
}


template < class T >
std::ostream& operator<< (std::ostream& st, const vec< T >& v)
{
	for (size_t i = 0; i < v.size(); i++)
	{
		st << v[i] << " ";
	}
	return st;
}


} /* end of namespace aml */

#endif /* _freeAML_vec_h_ */
