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


#ifndef _freeAML_mat_h_
#define _freeAML_mat_h_


#include <algorithm>

#include "vec.h"


namespace aml
{


template < class T >
class mat
{

public:

	typedef typename vec< T >::value_type           value_type;
	typedef typename vec< T >::reference            reference;
	typedef typename vec< T >::const_reference      const_reference;


	/***********************************************************************
	 *
	 *	CONSTRUCTORS
	 *
	 **********************************************************************/

	/** @brief default constructor */
	explicit mat (const size_t rows = 0,
	              const size_t cols = 0,
	              const T& defval = T());

	/** @brief copy constructor */
	mat (const mat< T >& A);


	/***********************************************************************
	 *
	 *	OVERLOADED OPERATORS
	 *
	 **********************************************************************/

	/** @brief gets a reference to element (i,j) of the matrix */
	reference operator() (const size_t i, const size_t j);

	/** @brief gets a const reference to element (i,j) of the matrix */
	const_reference operator() (const size_t i, const size_t j) const;

	/** @brief gets a reference to element (i,j) (passed as a 2d array) of the matrix */
	template < class _array2d >
	reference operator() (const _array2d& u);

	/** @brief gets a const reference to element (i,j) (passed as a 2d array) of the matrix */
	template < class _array2d >
	const_reference operator() (const _array2d& u) const;

	/** @brief gets a const reference to the i-th row of the matrix */
	const vec< T >& operator[] (const size_t i) const;

	/** @brief copies the elements of an input matrix */
	mat< T >& operator= (const mat< T >& B);

	/** @brief multiplies all elements of the matrix by a scalar */
	mat< T >& operator*= (const T& c);

	/** @brief computes the multiplication of the matrix by a scalar (on the right) */
	mat< T > operator* (const T& c) const;

	/** @brief divides all elements of the matrix by a scalar */
	mat< T >& operator/= (const T& c);

	/** @brief computes the division of the matrix by a scalar */
	mat< T > operator/ (const T& c) const;

	/** @brief adds an input matrix to this matrix */
	mat< T >& operator+= (const mat< T >& B);

	/** @brief computes the addition of this matrix and an input matrix */
	mat< T > operator+ (const mat< T >& B) const;

	/** @brief subtracts an input matrix from this matrix */
	mat< T >& operator-= (const mat< T >& B);

	/** @brief computes the subtraction of an input matrix from this matrix */
	mat< T > operator- (const mat< T >& B) const;

	/** @brief computes the multiplication of this matrix and an input matrix */
	mat< T > operator* (const mat< T >& B) const;

	/**
	 * @brief computes the multiplication of this matrix and a column vector
	 * @param v a vector (interpreted as a "column" vector)
	 * @return a vector which is the multiplication of this matrix and v
	 */
	vec< T > operator* (const vec< T >& v) const;


	/***********************************************************************
	 *
	 *	MATHEMATICAL FUNCTIONS
	 *
	 **********************************************************************/

	/** @brief computes the max norm of the matrix */
	T max_norm () const;

	/** @brief if this matrix is A, gets A(i_min:imax, j_min:j_max) */
	mat< T > submatrix(const size_t i_min,
	                   const size_t i_max,
	                   const size_t j_min,
	                   const size_t j_max) const;

	/** @brief computes the transpose of the matrix */
	mat< T > transpose () const;

	/** @brief returns true if the matrix is symmetric or false otherwise */
	bool is_symmetric () const;

	/** @brief returns true if the matrix is square or false otherwise */
	bool is_square () const;


	/***********************************************************************
	 *
	 *	OTHER FUNCTIONS
	 *
	 **********************************************************************/

	/** @brief sets all elements of the matrix to a constant */
	void fill (const T& x);

	/** @brief sets all elements of the matrix to zero */
	void zero_fill ();

	/** @brief gets the number of rows in the matrix */
	size_t num_rows () const;

	/** @brief gets the number of columns in the matrix */
	size_t num_cols () const;

	/** @brief returns true if the matrix is empty or false otherwise */
	bool empty () const;

	/** @brief resizes the matrix (new elements get an optional default value) */
	void resize (const size_t rows, const size_t cols, const T& defval = T());

	/** @brief swap rows i and j of the matrix */
	void swap_rows (const size_t i, const size_t j);


	/***********************************************************************
	 *
	 *	PRINT FUNCTIONS
	 *
	 **********************************************************************/

	/**
	 * @brief prints all elements of the matrix
	 * @param st the (optional) output stream
	 * @return the output stream st
	 */
	std::ostream& print (std::ostream& st = std::cout) const;

	/**
	 * @brief prints all elements of the matrix with their indices (i,j)
	 * @param st the (optional) output stream
	 * @return the output stream st
	 */
	std::ostream& print_with_indices (std::ostream& st = std::cout) const;


private:

	/* the matrix elements are stored as a vector containing its rows */
	vec< vec< T > > m_rows;

	/** @brief gets a reference to the i-th row of the matrix */
	vec< T >& operator[] (const size_t i);

}; /* end of class mat */


/*******************************************************************************
 *
 *	NON-MEMBER OVERLOADED OPERATORS AND FUNCTIONS
 *
 ******************************************************************************/

/**
 * @brief computes the multiplication of a "row" vector and a matrix
 * @param v a vector (interpreted as a "row" vector)
 * @param A a matrix
 * @return a vector which is the multiplication of v and A
 */
template < class T >
vec< T > operator* (const vec< T >& v, const mat< T >& A);

/** @brief computes the multiplication of a matrix by scalar (on the left) */
template < class T >
mat< T > operator* (const T& c, const mat< T >& A);

/** @brief prints the matrix to a given stream (st) */
template < class T >
std::ostream& operator<< (std::ostream& st, const mat< T >& A);

/** @brief returns the m × m identity matrix (stored as a dense matrix) */
template < class T >
mat< T > dense_identity_matrix (const size_t m);


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS AND MEMBER OPERATORS
 *
 ******************************************************************************/

template < class T >
mat< T >::mat (const size_t rows,
               const size_t cols,
               const T& defval): m_rows(rows, vec< T >(cols,defval))
{
	/* nothing needs to be done here */
}


template < class T >
mat< T >::mat (const mat< T >& A): m_rows(A.m_rows)
{
	/* nothing needs to be done here */
}


template < class T >
typename mat< T >::reference mat< T >::operator() (const size_t i, const size_t j)
{
	FREEAML_ASSERT(0 <= i &&  i < num_rows());
	FREEAML_ASSERT(0 <= j &&  j < num_cols());

	return m_rows[i][j];
}


template < class T >
typename mat< T >::const_reference mat< T >::operator() (const size_t i,
                                                         const size_t j) const
{
	FREEAML_ASSERT(0 <= i &&  i < num_rows());
	FREEAML_ASSERT(0 <= j &&  j < num_cols());

	return m_rows[i][j];
}


template < class T >
template < class _array2d >
typename mat< T >::reference mat< T >::operator() (const _array2d& u)
{
	return m_rows[u[0]][u[1]];
}


template < class T >
template < class _array2d >
typename mat< T >::const_reference mat< T >::operator() (const _array2d& u) const
{
	return m_rows[u[0]][u[1]];
}


template < class T >
mat< T >& mat< T >::operator= (const mat< T >& B)
{
	m_rows.resize(B.num_rows());

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < B.num_rows(); i++)
	{
		m_rows[i] = B[i];
	}

	return (*this);
}


template < class T >
const vec< T >& mat< T >::operator[] (const size_t i) const
{
	return m_rows[i];
}


template < class T >
mat< T >& mat< T >::operator*= (const T& c)
{

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < num_rows(); i++)
	{
		(*this)[i] *= c;
	}
	return (*this);
}


template < class T >
mat< T > mat< T >::operator* (const T& c) const
{
	mat< T > B(*this);
	B *= c;
	return B;
}


template < class T >
mat< T >& mat< T >::operator/= (const T& c)
{

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < num_rows(); i++)
	{
		(*this)[i] /= c;
	}
	return (*this);
}


template < class T >
mat< T > mat< T >::operator/ (const T& c) const
{
	mat< T > B(*this);
	B /= c;
	return B;
}


template < class T >
mat< T >& mat< T >::operator+= (const mat< T >& B)
{
	FREEAML_ASSERT(num_rows() == B.num_rows());
	FREEAML_ASSERT(num_cols() == B.num_cols());

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < num_rows(); i++)
	{
		(*this)[i] += B[i];
	}
	return (*this);
}


template < class T >
mat< T > mat< T >::operator+ (const mat< T >& B) const
{
	mat< T > C(*this);
	C += B;
	return C;
}


template < class T >
mat< T >& mat< T >::operator-= (const mat< T >& B)
{
	FREEAML_ASSERT(num_rows() == B.num_rows());
	FREEAML_ASSERT(num_cols() == B.num_cols());

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < num_rows(); i++)
	{
		(*this)[i] -= B[i];
	}
	return (*this);
}


template < class T >
mat< T > mat< T >::operator- (const mat< T >& B) const
{
	mat< T > C(*this);
	C -= B;
	return C;
}


template < class T >
mat< T > mat< T >::operator* (const mat< T >& B) const
{
	FREEAML_ASSERT(num_cols() == B.num_rows());

	mat< T > C(num_rows(), B.num_cols(), (T) 0);

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < num_rows(); i++)
	{
		for (size_t k = 0; k < num_cols(); k++)
		{
			for (size_t j = 0; j < B.num_cols(); j++)
			{
				/* C_ij = A_ik*B_kj (Einstein notation, A = this) */
				C(i,j) += (*this)(i,k) * B(k,j);
			}
		}
	}
	return C;
}


template < class T >
vec< T > mat< T >::operator* (const vec< T >& v) const
{
	FREEAML_ASSERT(num_cols() == v.size());

	vec< T > Av(num_rows(), (T) 0);

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < num_rows(); i++)
	{
		Av[i] = (*this)[i] * v;
	}
	return Av;
}


template < class T >
T mat< T >::max_norm () const
{
	T norm = (T) 0;

#ifdef _OPENMP
#pragma omp parallel
	{
		T local_norm = (T) 0;

#pragma omp for nowait

		for (size_t i = 0; i < num_rows(); i++)
		{
			/* this should also work with complex numbers */
			local_norm = std::max(std::abs(local_norm), (*this)[i].linf_norm());
		}

#pragma omp critical
		{
			norm = std::max(std::abs(local_norm), std::abs(norm));
		}
	}
#else

	/* serial implementation */
	for (size_t i = 0; i < num_rows(); i++)
	{
		norm = std::max(std::abs(norm), (*this)[i].linf_norm());
	}
#endif

	return norm;
}


template < class T >
mat< T > mat< T >::submatrix (const size_t i_min,
                              const size_t i_max,
                              const size_t j_min,
                              const size_t j_max) const
{
	FREEAML_ASSERT(i_min <= i_max && i_max < num_rows());
	FREEAML_ASSERT(j_min <= j_max && j_max < num_cols());

	mat< T > B(i_max-i_min+1, j_max-j_min+1);

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = i_min; i <= i_max; i++)
	{
		for (size_t j = j_min; j <= j_max; j++)
		{
			B(i-i_min, j-j_min) = (*this)(i,j);
		}
	}
	return B;
}


template < class T >
mat< T > mat< T >::transpose () const
{
	mat< T > At(num_cols(), num_rows(), T());

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < num_rows(); i++)
	{
		for (size_t j = 0; j < num_cols(); j++)
		{
			At(j,i) = (*this)(i,j);
		}
	}
	return At;
}


template < class T >
bool mat< T >::is_symmetric () const
{
	if (!is_square())
	{
		return false;
	}

	for (size_t i = 0; i < num_rows(); i++)
	{
		for (size_t j = i+1; j < num_cols(); j++)
		{
			if ((*this)(i,j) != (*this)(j,i))
			{
				return false;
			}
		}
	}
	return true;
}


template < class T >
bool mat< T >::is_square () const
{
	return (num_rows() == num_cols());
}


template < class T >
void mat<T>::fill(const T& x)
{

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < num_rows(); i++)
	{
		(*this)[i].fill(x);
	}
}


template < class T >
void mat<T>::zero_fill()
{
	fill((T)0);
}


template < class T >
size_t mat< T >::num_rows () const
{
	return m_rows.size();
}


template < class T >
size_t mat< T >::num_cols () const
{
	return m_rows.empty() ? 0 : m_rows[0].size();
}


template < class T >
bool mat< T >::empty () const
{
	return m_rows.empty();
}


template < class T >
void mat<T>::resize(const size_t rows, const size_t cols, const T& defval)
{
	FREEAML_ASSERT(rows >= 0);
	FREEAML_ASSERT(cols >= 0);

	m_rows.resize(rows);

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < rows; i++)
	{
		m_rows[i].resize(cols,defval);
	}
}


template < class T >
void mat< T >::swap_rows (const size_t i, const size_t j)
{
	FREEAML_ASSERT(0 <= i && i < num_rows());
	FREEAML_ASSERT(0 <= j && j < num_rows());

	std::swap(m_rows[i],m_rows[j]);
}


template < class T >
std::ostream& mat< T >::print (std::ostream& st) const
{
	size_t max_width = 0;

	for (size_t i = 0; i < num_rows(); i++)
	{
		for (size_t j = 0; j < num_cols(); j++)
		{
			std::stringstream ss;
			ss << (*this)(i,j);

			max_width = std::max(max_width, ss.str().size());
		}
	}

	for (size_t i = 0; i < num_rows(); i++)
	{
		for (size_t j = 0; j < num_cols(); j++)
		{
			st << std::setw(max_width+2) << (*this)(i,j);
		}
		st << "\n";
	}

	return st;
}


template < class T >
std::ostream& mat< T >::print_with_indices (std::ostream& st) const
{
	size_t max_width = 0;

	for (size_t i = 0; i < num_rows(); i++)
	{
		for (size_t j = 0; j < num_cols(); j++)
		{
			std::stringstream ss;
			ss << (*this)(i,j);

			max_width = std::max(max_width, ss.str().size());
		}
	}

	size_t nr = num_digits(num_rows());
	size_t nc = num_digits(num_cols());

	for (size_t i = 0; i < num_rows(); i++)
	{
		for (size_t j = 0; j < num_cols(); j++)
		{
			st << "(" << std::setw(nr) << i
			   << "," << std::setw(nc) << j
			   << ")" << " "
			   << std::setw(max_width+2)
			   << (*this)(i,j)
			   << "  ";
		}
		st << "\n";
	}

	return st;
}


template < class T >
vec< T >& mat< T >::operator[] (const size_t i)
{
	return m_rows[i];
}


/*******************************************************************************
 *
 *	DEFINITIONS OF NON-MEMBER OVERLOADED OPERATORS AND FUNCTIONS
 *
 ******************************************************************************/

template < class T >
vec< T > operator* (const vec< T >& v, const mat< T >& A)
{
	FREEAML_ASSERT(v.size() == A.num_rows());

	vec< T > vA(A.num_cols(), (T) 0);

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t j = 0; j < A.num_cols(); j++)
	{
		for (size_t i = 0; i < A.num_rows(); i++)
		{
			vA[j] += v[i] * A(i,j);
		}
	}
	return vA;
}


template < class T >
mat< T > operator* (const T& c, const mat< T >& A)
{
	return A*c;
}


template < class T >
std::ostream& operator<< (std::ostream& st, const mat< T >& A)
{
	return A.print(st);
}


template < class T >
mat< T > dense_identity_matrix (const size_t m)
{
	mat< T > A(m, m, (T) 0);

	for (size_t i = 0; i < m; i++)
	{
		A(i,i) = (T) 1;
	}
	return A;
}


} /* end of namespace aml */

#endif /* _freeAML_mat_h_ */
