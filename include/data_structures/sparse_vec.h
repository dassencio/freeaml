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


#ifndef _freeAML_sparse_vec_h_
#define _freeAML_sparse_vec_h_


#include <map>

#include "vec.h"


namespace aml
{


template < class T >
class sparse_vec
{

public:

	typedef std::map< size_t, T >                   map_type;
	typedef std::map< const size_t, T >             const_map_type;
	typedef typename map_type::value_type           value_type;
	typedef typename map_type::iterator             iterator;
	typedef typename map_type::const_iterator       const_iterator;


	/***********************************************************************
	 *
	 *	CONSTRUCTORS
	 *
	 **********************************************************************/

	/** @brief default constructor */
	explicit sparse_vec (const size_t len = 0);

	/** @brief range constructor */
	template < class _iterator >
	sparse_vec (_iterator first, _iterator last);

	/** @brief copy constructor */
	sparse_vec (const sparse_vec< T >& v);


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
	sparse_vec< T >& operator*= (const T& c);

	/** @brief computes the multiplication of the vector by a scalar (on the right) */
	sparse_vec< T > operator* (const T& c) const;

	/** @brief divides all elements of the vector by a scalar */
	sparse_vec< T >& operator/= (const T& c);

	/** @brief computes the division of the vector by a scalar */
	sparse_vec< T > operator/ (const T& c) const;

	/** @brief adds an input vector to this vector */
	sparse_vec< T >& operator+= (const sparse_vec< T >& v);

	/** @brief computes the addition of this vector and an input vector */
	sparse_vec< T > operator+ (const sparse_vec< T >& v) const;

	/** @brief subtracts an input vector from this vector */
	sparse_vec< T >& operator-= (const sparse_vec< T >& v);

	/** @brief computes the subtraction of an input vector from this vector */
	sparse_vec< T > operator- (const sparse_vec< T >& v) const;

	/** @brief computes the dot product of this vector and an input sparse vector */
	T operator* (const sparse_vec< T >& v) const;

	/** @brief computes the dot product of this vector and an input dense vector */
	T operator* (const vec< T >& v) const;

	/** @brief returns true if this vector is equal to an input vector */
	bool operator== (const sparse_vec< T >& v) const;

	/** @brief returns true if this vector is not equal to an input vector */
	bool operator!= (const sparse_vec< T >& v) const;


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


	/***********************************************************************
	 *
	 *	OTHER FUNCTIONS
	 *
	 **********************************************************************/

	/** @brief sets all elements of the vector to zero */
	sparse_vec< T >& zero_fill ();

	/** @brief computes the number of nonzero elements stored in the vector */
	size_t num_nonzero_elements () const;

	/** @brief gets the total number of elements stored in the vector (if
	 *         zeros were stored by hand, they will also be counted!) */
	size_t num_stored_elements () const;

	/**
	 * @brief removes all elements of the vector storing zeros
	 * @note the vector size remains the same; this function merely
	 *       frees the memory associated with elements storing zeros
	 *       (however, avoid adding zeros to the sparse vector by only
	 *       reading its elements in const mode and by not writing zeros to
	 *       it in the first place!)
	 */
	void remove_zeros();

	/** @brief gets the size of the vector */
	size_t size () const;

	/** @brief returns true if the vector is empty or false otherwise */
	bool empty () const;

	/** @brief clears the vector */
	void clear ();

	/** @brief resizes the vector (new elements are assumed to be zero) */
	void resize (const size_t len);

	/** @brief inserts an element at the end of the vector */
	void push_back (const T& x);

	/** @brief removes the last element of the vector */
	void pop_back ();

	/** @brief gets an iterator to the beginning of the vector */
	iterator begin ();

	/** @brief gets a const iterator to the beginning of the vector */
	const_iterator begin() const;

	/** @brief gets an iterator to the end of the vector */
	iterator end ();

	/** @brief gets a const iterator to the end of the vector */
	const_iterator end () const;


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


private:

	/* a map stores all elements of the sparse vector */
	map_type	m_elements;

	/* length (size) of the sparse vector (this value can be differnt
	 * from the number of elements actually stored in the sparse vector) */
	size_t		m_len;

}; /* end of class sparse_vec */


/*******************************************************************************
 *
 *	NON-MEMBER OVERLOADED OPERATORS
 *
 ******************************************************************************/

/** @brief computes the multiplication of a sparse vector by a scalar */
template < class T >
sparse_vec< T > operator* (const T& c, const sparse_vec< T >& v);

/** @brief computes the dot product of a dense vector and a sparse vector */
template < class T >
T operator* (const vec< T >& u, const sparse_vec< T >& v);

/** @brief prints the elements of a sparse vector directly to an output stream */
template < class T >
std::ostream& operator<< (std::ostream& st, const sparse_vec< T >& v);


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS AND MEMBER OPERATORS
 *
 ******************************************************************************/

template < class T >
sparse_vec< T >::sparse_vec (const size_t len): m_len(len)
{
	/* nothing needs to be done here */
}


template < class T >
template < class _iterator >
sparse_vec< T >::sparse_vec (_iterator first, _iterator last): m_elements(first,last)
{
	/* nothing needs to be done here */
}


template < class T >
sparse_vec< T >::sparse_vec (const sparse_vec< T >& v):

	m_elements(v.m_elements), m_len(v.m_len)
{
	/* nothing needs to be done here */
}


template < class T >
T& sparse_vec< T >::operator[] (const size_t i)
{
	FREEAML_ASSERT(i < size());

	return m_elements[i];
}


template < class T >
const T& sparse_vec< T >::operator[] (const size_t i) const
{
	FREEAML_ASSERT(i < size());

	static T m_zero = (T) 0;

	const_iterator it;

	if ((it = m_elements.find(i)) != end())
	{
		return it->second;
	}
	return m_zero;
}


template < class T >
sparse_vec< T >& sparse_vec< T >::operator*= (const T& c)
{
	for (iterator it = begin(); it != end(); it++)
	{
		it->second *= c;
	}
	return (*this);
}


template < class T >
sparse_vec< T > sparse_vec< T >::operator* (const T& c) const
{
	sparse_vec< T > u = (*this);
	return (u *= c);
}


template < class T >
sparse_vec< T >& sparse_vec< T >::operator/= (const T& c)
{
	for (iterator it = begin(); it != end(); it++)
	{
		it->second /= c;
	}
	return (*this);
}


template < class T >
sparse_vec< T > sparse_vec< T >::operator/ (const T& c) const
{
	sparse_vec< T > u = (*this);
	return (u /= c);
}


template < class T >
sparse_vec< T >& sparse_vec< T >::operator+= (const sparse_vec< T >& v)
{
	FREEAML_ASSERT(size() == v.size());

	for (const_iterator it = v.begin(); it != v.end(); it++)
	{
		m_elements[it->first] += it->second;
	}
	return (*this);
}


template < class T >
sparse_vec< T > sparse_vec< T >::operator+ (const sparse_vec< T >& v) const
{
	FREEAML_ASSERT(size() == v.size());

	sparse_vec< T > u = (*this);
	return (u += v);
}


template < class T >
sparse_vec< T >& sparse_vec< T >::operator-= (const sparse_vec< T >& v)
{
	FREEAML_ASSERT(size() == v.size());

	for (const_iterator it = v.begin(); it != v.end(); it++)
	{
		m_elements[it->first] -= it->second;
	}
	return (*this);
}


template < class T >
sparse_vec< T > sparse_vec< T >::operator- (const sparse_vec< T >& v) const
{
	FREEAML_ASSERT(size() == v.size());

	sparse_vec< T > u = (*this);
	return (u -= v);
}


template < class T >
T sparse_vec< T >::operator* (const sparse_vec< T >& v) const
{
	FREEAML_ASSERT(size() == v.size());

	T dotp = (T) 0;

	for (const_iterator it = begin(); it != end(); it++)
	{
		dotp += it->second * v[it->first];
	}
	return dotp;
}


template < class T >
T sparse_vec< T >::operator* (const vec< T >& v) const
{
	FREEAML_ASSERT(size() == v.size());

	T dotp = (T) 0;

	for (const_iterator it = begin(); it != end(); it++)
	{
		dotp += it->second * v[it->first];
	}
	return dotp;
}


template < class T >
bool sparse_vec< T >::operator== (const sparse_vec< T >& v) const
{
	return (m_elements == v.m_elements && size() == v.size());
}


template < class T >
bool sparse_vec< T >::operator!= (const sparse_vec< T >& v) const
{
	return !(*this == v);
}


template < class T >
T sparse_vec< T >::l1_norm () const
{
	T norm = (T) 0;

	for (const_iterator it = begin(); it != end(); it++)
	{
		norm += std::abs(it->second);
	}
	return norm;
}


template < class T >
T sparse_vec< T >::l2_norm () const
{
	T norm = (T) 0;

	for (const_iterator it = begin(); it != end(); it++)
	{
		norm += std::abs(it->second) * std::abs(it->second);
	}
	return std::sqrt(norm);
}


template < class T >
T sparse_vec< T >::lp_norm (const int p) const
{
	FREEAML_ASSERT(p > 0);

	T norm = (T) 0;

	for (const_iterator it = begin(); it != end(); it++)
	{
		norm += std::pow(std::abs(it->second), p);
	}
	return std::pow(norm, (T)1/p);
}


template < class T >
T sparse_vec< T >::linf_norm () const
{
	T norm = (T) 0;

	for (const_iterator it = begin(); it != end(); it++)
	{
		/* this should also work with complex numbers */
		norm = std::max(std::abs(norm), std::abs(it->second));
	}
	return norm;
}


template < class T >
T sparse_vec< T >::sum (void) const
{
	T sum = (T) 0;

	for (const_iterator it = begin(); it != end(); it++)
	{
		sum += it->second;
	}
	return sum;
}


template < class T >
T sparse_vec< T >::average () const
{
	if (!empty())
	{
		return sum() / (T) size();
	}
	return (T) 0;
}


template < class T >
sparse_vec< T >& sparse_vec< T >::zero_fill ()
{
	m_elements.clear();
	return (*this);
}


template < class T >
size_t sparse_vec< T >::num_nonzero_elements () const
{
	size_t count = 0;

	for (const_iterator it = begin(); it != end(); it++)
	{
		if (it->second != (T) 0)
		{
			count++;
		}
	}
	return count;
}


template < class T >
size_t sparse_vec< T >::num_stored_elements () const
{
	return m_elements.size();
}


template < class T >
void sparse_vec< T >::remove_zeros()
{
	iterator it = begin();

	while (it != end())
	{
		if (it->second == (T) 0)
		{
			/* safe erase */
			m_elements.erase(it++);
		}
		else
		{
			it++;
		}
	}
}


template < class T >
size_t sparse_vec< T >::size () const
{
	return m_len;
}


template < class T >
bool sparse_vec< T >::empty () const
{
	return (size() == 0);
}


template < class T >
void sparse_vec<T>::clear()
{
	m_elements.clear();
	m_len = 0;
	return (*this);
}


template < class T >
void sparse_vec< T >::resize (const size_t len)
{
	FREEAML_ASSERT(len >= 0);

	if (len < size())
	{
		/* start at the end of the sparse vector (largest key) */
		iterator it = end();

		while (it != begin())
		{
			if (it->first >= len)
			{
				/* safe erase */
				m_elements.erase(it--);
			}
			else
			{
				/* the vector elements are ordered by key */
				break;
			}
		}
	}
	m_len = len;
}


template < class T >
void sparse_vec< T >::push_back (const T& x)
{
	m_elements[size()] = x;
	m_len++;
	return (*this);
}


template < class T >
void sparse_vec< T >::pop_back ()
{
	if (!empty())
	{
		m_elements.erase(end());
		m_len--;
	}
	return (*this);
}


template < class T >
typename sparse_vec< T >::iterator sparse_vec< T >::begin ()
{
	return m_elements.begin();
}


template < class T >
typename sparse_vec< T >::const_iterator sparse_vec< T >::begin () const
{
	return m_elements.begin();
}



template < class T >
typename sparse_vec< T >::iterator sparse_vec< T >::end ()
{
	return m_elements.end();
}


template < class T >
typename sparse_vec< T >::const_iterator sparse_vec< T >::end () const
{
	return m_elements.end();
}


template < class T >
std::ostream& sparse_vec< T >::print (std::ostream& st) const
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
std::ostream& sparse_vec< T >::print_with_indices (std::ostream& st) const
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
sparse_vec< T > operator* (const T& c, const sparse_vec< T >& v)
{
	return (v*c);
}


template < class T >
T operator* (const vec< T >& u, const sparse_vec< T >& v)
{
	return (v*u);
}


template < class T >
std::ostream& operator<< (std::ostream& st, const sparse_vec< T >& v)
{
	for (size_t i = 0; i < v.size(); i++)
	{
		st << v[i] << " ";
	}
	return st;
}


} /* end of namespace aml */

#endif /* _freeAML_sparse_vec_h_ */
