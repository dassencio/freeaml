#ifndef _freeAML_generic_iter_h_
#define _freeAML_generic_iter_h_


#include <cstddef>

#include <general/debug.h>


/*******************************************************************************
 *
 *	GENERIC ITERATIVE LINEAR SYSTEM SOLVER (ABSTRACT CLASS)
 *
 ******************************************************************************/


namespace aml
{


template < class T >
class generic_iter
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
	generic_iter (const size_t max_iter, const T& residual_tol);


	/***********************************************************************
	 *
	 *	SOLVE FUNCTIONS
	 *
	 **********************************************************************/

	/* a solve function must be implemented in each iterative method */


	/***********************************************************************
	 *
	 *	CONFIGURATION FUNCTIONS
	 *
	 **********************************************************************/

	/** @brief gets the number of iterations performed */
	size_t num_iterations () const;

	/** @brief gets the maximum number of iterations allowed in one solve */
	size_t max_num_iterations () const;

	/** @brief gets the maximum residual tolerance */
	T residual_tolerance () const;

	/** @brief sets the number of iterations performed */
	void set_num_iterations (const size_t num_iter);

	/** @brief sets the maximum number of iterations allowed in one solve */
	void set_max_num_iterations (const size_t max_iter);

	/** @brief sets the maximum residual tolerance */
	void set_residual_tolerance (const T& residual_tol);


	/***********************************************************************
	 *
	 *	OTHER FUNCTIONS
	 *
	 **********************************************************************/

	/** @brief checks if the dimensions of the system matrix, the solution
	 *         vector and the right-hand side vector match */
	template < class _mat, class _vec >
	bool check_dimensions (const _mat& A, const _vec& x, const _vec& b) const;


private:

	size_t          m_num_iter;
	size_t          m_max_iter;
	T               m_residual_tol;

}; /* end of class generic_iter */


/*******************************************************************************
 *
 *	DEFINITIONS OF CLASS MEMBER FUNCTIONS
 *
 ******************************************************************************/

template < class T >
generic_iter< T >::generic_iter (const size_t max_iter, const T &residual_tol):

        m_num_iter(0), m_max_iter(max_iter), m_residual_tol(residual_tol)
{
	FREEAML_ASSERT(max_iter > 0);
	FREEAML_ASSERT(residual_tol > (T) 0);

	/* nothing needs to be done here */
}


template < class T >
size_t generic_iter< T >::num_iterations () const
{
	return m_num_iter;
}


template < class T >
size_t generic_iter< T >::max_num_iterations () const
{
	return m_max_iter;
}


template < class T >
T generic_iter< T >::residual_tolerance () const
{
	return m_residual_tol;
}


template < class T >
void generic_iter< T >::set_num_iterations (const size_t num_iter)
{
	m_num_iter = num_iter;
}


template < class T >
void generic_iter< T >::set_max_num_iterations (const size_t max_iter)
{
	m_max_iter =  max_iter;
}


template < class T >
void generic_iter< T >::set_residual_tolerance (const T& residual_tol)
{
	FREEAML_ASSERT(residual_tol > (T) 0);

	m_residual_tol = residual_tol;
}


template < class T >
template < class _mat, class _vec >
bool generic_iter< T >::check_dimensions (const _mat& A,
                                          const _vec& x,
                                          const _vec& b) const
{
	return (A.num_rows() > 0) &&
	               A.is_square() &&
	                   x.size() == b.size() &&
	                       A.num_cols() == x.size();
}


} /* end of namespace aml */


#endif /* _freeAML_generic_iter_h_ */
