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


#ifndef _freeAML_sparse_mat_h_
#define _freeAML_sparse_mat_h_


#include "mat.h"
#include "sparse_vec.h"


namespace aml
{


template < class T >
class sparse_mat
{

public:

	typedef typename sparse_vec< T >::value_type		value_type;

	/* the scope of these iterators is a single row, not the entire
	 * matrix (yes, the names are somewhat abused in this context) */
	typedef typename sparse_vec< T >::iterator		iterator;
	typedef typename sparse_vec< T >::const_iterator	const_iterator;


	/***********************************************************************
	 *
	 *	CONSTRUCTORS
	 *
	 **********************************************************************/

	/** @brief default constructor */
	explicit sparse_mat (const size_t rows = 0, const size_t cols = 0);

	/** @brief copy constructor */
	sparse_mat (const sparse_mat< T >& A);


	/***********************************************************************
	 *
	 *	OVERLOADED OPERATORS
	 *
	 **********************************************************************/

	/** @brief gets a reference to element (i,j) of the matrix */
	T& operator() (const size_t i, const size_t j);

	/** @brief gets a const reference to element (i,j) of the matrix */
	const T& operator() (const size_t i, const size_t j) const;

	/** @brief gets a reference to element (i,j) (passed as a 2d array) of the matrix */
	template < class _array2d >
	T& operator() (const _array2d& u);

	/** @brief gets a const reference to element (i,j) (passed as a 2d array) of the matrix */
	template < class _array2d >
	const T& operator() (const _array2d& u) const;

	/** @brief gets a const reference to the i-th row of the matrix */
	const sparse_vec< T >& operator[] (const size_t i) const;

	/** @brief multiplies all elements of the matrix by a scalar */
	sparse_mat< T >& operator*= (const T& c);

	/** @brief computes the multiplication of the matrix by a scalar (on the right) */
	sparse_mat< T > operator* (const T& c) const;

	/** @brief divides all elements of the matrix by a scalar */
	sparse_mat< T >& operator/= (const T& c);

	/** @brief computes the division of the matrix by a scalar */
	sparse_mat< T > operator/ (const T& c) const;

	/** @brief adds an input matrix to this matrix */
	sparse_mat< T >& operator+= (const sparse_mat< T >& B);

	/** @brief computes the addition of this matrix and an input matrix */
	sparse_mat< T > operator+ (const sparse_mat< T >& B) const;

	/** @brief subtracts an input matrix from this matrix */
	sparse_mat< T >& operator-= (const sparse_mat< T >& B);

	/** @brief computes the subtraction of an input matrix from this matrix */
	sparse_mat< T > operator- (const sparse_mat< T >& B) const;

	/** @brief computes the multiplication of this matrix and an input sparse matrix */
	sparse_mat< T > operator* (const sparse_mat< T >& B) const;

	/**
	 * @brief computes the multiplication of this matrix and an input dense matrix
	 * @return the matrix product as a dense matrix
	 */
	mat< T > operator* (const mat< T >& B) const;

	/**
	 * @brief computes the multiplication of this matrix and a dense "column" vector
	 * @param v a dense vector (interpreted as a "column" vector)
	 * @return a dense vector which is the multiplication of this matrix and v
	 */
	vec< T > operator* (const vec< T >& v) const;


	/***********************************************************************
	 *
	 *	MATHEMATICAL FUNCTIONS
	 *
	 **********************************************************************/

	/** @brief computes the max norm of the matrix */
	T max_norm () const;

	/** @brief computes the transpose of the matrix */
	sparse_mat< T > transpose () const;

	/** @brief returns true if the matrix is symmetric or false otherwise */
	bool is_symmetric () const;

	/** @brief returns true if the matrix is square or false otherwise */
	bool is_square () const;


	/***********************************************************************
	 *
	 *	OTHER FUNCTIONS
	 *
	 **********************************************************************/

	/** @brief sets all values of the matrix to zero */
	void set_all_values_to_zero ();

	/** @brief computes the number of nonzero values stored on the matrix */
	size_t num_nonzero_values () const;

	/** @brief gets the total number of values stored on the matrix (if
	 *	   zeros were stored by hand, they will also be counted!) */
	size_t num_stored_values () const;

	/**
	 * @brief removes all elements of the matrix storing zeros
	 * @note the matrix dimensions remain the same; this function merely
	 *       frees the memory associated with elements storing zeros
	 *       (however, avoid adding zeros to the sparse matrix by only
	 *       reading its values in const mode and by not writing zeros to
	 *       it in the first place!)
	 */
	void remove_zeros();

	/** @brief gets the number of rows in the matrix */
	size_t num_rows () const;

	/** @brief gets the number of columns in the matrix */
	size_t num_cols () const;

	/** @brief returns true if the matrix is empty or false otherwise */
	bool empty () const;

	/** @brief clears the matrix */
	void clear ();

	/** @brief resizes the matrix */
	void resize (const size_t rows, const size_t cols);

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

	/**
	 * @brief prints all nonzero elements of the matrix with their indices (i,j)
	 * @param width width of each printed element
	 * @param prec precision of each printed element (for doubles, floats...)
	 * @param st output stream
	 * @return the output stream st
	 */
	std::ostream& print_nonzero_with_indices (std::ostream& st = std::cout) const;

private:

	/* the matrix elements are stored as a vector containing its rows */
	vec< sparse_vec< T > >	m_rows;

	/** @brief gets a reference to the i-th row of the matrix */
	sparse_vec< T >& operator[] (const size_t i);

}; /* end of class sparse_mat */


/*******************************************************************************
 *
 *	NON-MEMBER OVERLOADED OPERATORS AND FUNCTIONS
 *
 ******************************************************************************/

/**
 * @brief computes the multiplication of a dense matrix and a sparse matrix
 * @param A a dense matrix
 * @param B a sparse matrix
 * @return the matrix product of A and B as a dense matrix
 */
template < class T >
mat< T > operator* (const mat< T >& A, const sparse_mat< T >& B);

/**
 * @brief computes the multiplication of a dense "row" vector and a sparse matrix
 * @param v a vector (interpreted as a "row" vector)
 * @param A a sparse matrix
 * @return a dense vector which is the multiplication of v and A
 */
template < class T >
vec< T > operator* (const vec< T >& v, const sparse_mat< T >& A);

/** @brief computes the multiplication of a matrix by scalar (on the left) */
template < class T >
sparse_mat< T > operator* (const T& c, const sparse_mat< T >& A);

/** @brief returns the m × m identity matrix */
template < class T >
sparse_mat< T > identity_matrix (const size_t m);


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS AND MEMBER OPERATORS
 *
 ******************************************************************************/

template < class T >
sparse_mat< T >::sparse_mat (const size_t rows, const size_t cols):

	m_rows(rows, sparse_vec< T >(cols))
{
	/* nothing needs to be done here */
}


template < class T >
sparse_mat< T >::sparse_mat (const sparse_mat< T >& A): m_rows(A.m_rows)
{
	/* nothing needs to be done here */
}


template < class T >
T& sparse_mat< T >::operator() (const size_t i, const size_t j)
{
	FREEAML_ASSERT(0 <= i &&  i < num_rows());
	FREEAML_ASSERT(0 <= j &&  j < num_cols());

	return m_rows[i][j];
}


template < class T >
const T& sparse_mat< T >::operator() (const size_t i, const size_t j) const
{
	FREEAML_ASSERT(0 <= i &&  i < num_rows());
	FREEAML_ASSERT(0 <= j &&  j < num_cols());

	return m_rows[i][j];
}


template < class T >
template < class _array2d >
T& sparse_mat< T >::operator() (const _array2d& u)
{
	return m_rows[u[0]][u[1]];
}


template < class T >
template < class _array2d >
const T& sparse_mat< T >::operator() (const _array2d& u) const
{
	return m_rows[u[0]][u[1]];
}


template < class T >
const sparse_vec< T >& sparse_mat< T >::operator[] (const size_t i) const
{
	return m_rows[i];
}


template < class T >
sparse_mat< T >& sparse_mat< T >::operator*= (const T& c)
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
sparse_mat< T > sparse_mat< T >::operator* (const T& c) const
{
	sparse_mat< T > B(*this);
	B *= c;
	return B;
}


template < class T >
sparse_mat< T >& sparse_mat< T >::operator/= (const T& c)
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
sparse_mat< T > sparse_mat< T >::operator/ (const T& c) const
{
	sparse_mat< T > B(*this);
	B /= c;
	return B;
}


template < class T >
sparse_mat< T >& sparse_mat< T >::operator+= (const sparse_mat< T >& B)
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
sparse_mat< T > sparse_mat< T >::operator+ (const sparse_mat< T >& B) const
{
	sparse_mat< T > C(*this);
	C += B;
	return C;
}


template < class T >
sparse_mat< T >& sparse_mat< T >::operator-= (const sparse_mat< T >& B)
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
sparse_mat< T > sparse_mat< T >::operator- (const sparse_mat< T >& B) const
{
	sparse_mat< T > C(*this);
	C -= B;
	return C;
}


template < class T >
sparse_mat< T > sparse_mat< T >::operator* (const sparse_mat< T >& B) const
{
	FREEAML_ASSERT(num_cols() == B.num_rows());

	sparse_mat< T > C(num_rows(), B.num_cols());

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < num_rows(); i++)
	{
		for (const_iterator k = (*this)[i].begin(); k != (*this)[i].end(); k++)
		{
			for (const_iterator j = B[k->first].begin(); j != B[k->first].end(); j++)
			{
				/* C_ij = A_ik*B_kj (Einstein notation, A = this) */
				C(i,j->first) += k->second * j->second;
			}
		}
	}
	return C;
}


template < class T >
mat< T > sparse_mat< T >::operator* (const mat< T >& B) const
{
	FREEAML_ASSERT(num_cols() == B.num_rows());

	mat< T > C(num_rows(), B.num_cols(), (T)0);

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < num_rows(); i++)
	{
		for (const_iterator k = (*this)[i].begin(); k != (*this)[i].end(); k++)
		{
			for (size_t j = 0; j < B.num_cols(); j++)
			{
				/* C_ij = A_ik*B_kj (Einstein notation, A = this) */
				C(i,j) += k->second * B(k->first,j);
			}
		}
	}
	return C;
}


template < class T >
vec< T > sparse_mat< T >::operator* (const vec< T >& v) const
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
T sparse_mat< T >::max_norm () const
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
sparse_mat< T > sparse_mat< T >::transpose () const
{
	sparse_mat< T > At(num_cols(), num_rows());

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < num_rows(); i++)
	{
		for (const_iterator j = (*this)[i].begin(); j != (*this)[i].end(); j++)
		{
			At(j->first,i) = j->second;
		}
	}
	return At;
}


template < class T >
bool sparse_mat< T >::is_symmetric () const
{
	if (!is_square())
	{
		return false;
	}

	for (size_t i = 0; i < num_rows(); i++)
	{
		for (const_iterator j = (*this)[i].begin(); j != (*this)[i].end(); j++)
		{
			if (j->second != (*this)(j->first,i))
			{
				return false;
			}
		}
	}
	return true;
}


template < class T >
bool sparse_mat< T >::is_square () const
{
	return (num_rows() == num_cols());
}


template < class T >
void sparse_mat< T >::set_all_values_to_zero()
{

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < num_rows(); i++)
	{
		(*this)[i].set_all_values_to_zero();
	}
}


template < class T >
size_t sparse_mat< T >::num_nonzero_values () const
{
	size_t count = 0;

	for (size_t i = 0; i < num_rows(); i++)
	{
		count += (*this)[i].num_nonzero_values();
	}
	return count;
}


template < class T >
size_t sparse_mat< T >::num_stored_values () const
{
	size_t count = 0;

	for (size_t i = 0; i < num_rows(); i++)
	{
		count += (*this)[i].num_stored_values();
	}
	return count;
}


template < class T >
void sparse_mat< T >::remove_zeros()
{
	for (size_t i = 0; i < num_rows(); i++)
	{
		(*this)[i].remove_zeros();
	}
}


template < class T >
size_t sparse_mat< T >::num_rows () const
{
	return m_rows.size();
}


template < class T >
size_t sparse_mat< T >::num_cols () const
{
	return m_rows.empty() ? 0 : m_rows[0].size();
}


template < class T >
bool sparse_mat< T >::empty () const
{
	return m_rows.empty();
}


template < class T >
void sparse_mat< T >::clear ()
{
	m_rows.clear();
}


template < class T >
void sparse_mat<T>::resize(const size_t rows, const size_t cols)
{
	FREEAML_ASSERT(rows >= 0);
	FREEAML_ASSERT(cols >= 0);

	m_rows.resize(rows);

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < rows; i++)
	{
		m_rows[i].resize(cols);
	}
}


template < class T >
void sparse_mat< T >::swap_rows (const size_t i, const size_t j)
{
	FREEAML_ASSERT(0 <= i && i < num_rows());
	FREEAML_ASSERT(0 <= j && j < num_rows());

	std::swap(m_rows[i],m_rows[j]);
}


template < class T >
std::ostream& sparse_mat< T >::print (std::ostream& st) const
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
std::ostream& sparse_mat< T >::print_with_indices (std::ostream& st) const
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
std::ostream& sparse_mat< T >::print_nonzero_with_indices (std::ostream& st) const
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
			if ((*this)(i,j) != (T) 0)
			{
				st << "(" << std::setw(nr) << i
				   << "," << std::setw(nc) << j
				   << ")" << " "
				   << std::setw(max_width+2)
				   << (*this)(i,j)
				   << "  ";
			}
		}
		st << "\n";
	}

	return st;
}


template < class T >
sparse_vec< T >& sparse_mat< T >::operator[] (const size_t i)
{
	return m_rows[i];
}


/*******************************************************************************
 *
 *	DEFINITIONS OF NON-MEMBER OVERLOADED OPERATORS AND FUNCTIONS
 *
 ******************************************************************************/

template < class T >
sparse_mat< T > operator* (const mat< T >& A, const sparse_mat< T >& B)
{
	FREEAML_ASSERT(A.num_cols() == B.num_rows());

	mat< T > C(A.num_rows(), B.num_cols(), (T)0);

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < A.num_rows(); i++)
	{
		for (size_t k = 0; k < A.num_cols(); k++)
		{
			typedef typename sparse_mat< T >::const_iterator const_iterator;

			for (const_iterator j = B[k].begin(); j != B[k].end(); j++)
			{
				/* C_ij = A_ik*B_kj (Einstein notation) */
				C(i,j->first) += A(i,k) * j->second;
			}
		}
	}
	return C;
}


template < class T >
vec< T > operator* (const vec< T >& v, const sparse_mat< T >& A)
{
	FREEAML_ASSERT(v.size() == A.num_rows());

	vec< T > vA(A.num_cols(), (T) 0);

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (size_t i = 0; i < A.num_rows(); i++)
	{
		typedef typename sparse_mat< T >::const_iterator const_iterator;

		for (const_iterator j = A[i].begin(); j != A[i].end(); j++)
		{
			vA[j->first] += v[i] * j->second;
		}
	}
	return vA;
}


template < class T >
sparse_mat< T > operator* (const T& c, const sparse_mat< T >& A)
{
	return A*c;
}


template < class T >
sparse_mat< T > identity_matrix (const size_t m)
{
	sparse_mat< T > A(m, m);

	for (size_t i = 0; i < m; i++)
	{
		A(i,i) = (T) 1;
	}
	return A;
}


} /* end of namespace aml */

#endif /* _freeAML_sparse_mat_h_ */
