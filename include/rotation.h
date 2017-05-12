#pragma once

#include <debug.h>
#include <cmath>

namespace freeaml
{
/**
 * @brief Computes the Householder vector and the associated parameters which
 *        define the Householder transformation for a given input vector.
 * @param u The vector for which the Householder transformation will be defined.
 * @param v At the end, the Householder transformation vector for @c u.
 * @param b At the end, the first Householder transformation constant.
 * @param c At the end, the second Householder transformation constant.
 * @param j The index such that the Householder transformation for @c u projects
 *        it onto the unit vector <tt>e<sub>j</sub></tt> whose elements are all
 *        zero except for the <tt>j</tt>-th element, which is @c 1.
 * @note The vector @c v and the constants @c b and @c c are such that
 *       <tt>(I - bvv<sup>t</sup>)u = -c|u|e<sub>j</sub></tt>, where @c I is the
 *       identity matrix, @c v is a column vector and @c |u| is the
 *       @a L<sup>2</sup> norm of @c u.
 */
template<typename T, typename VectorType>
void householder(const VectorType& u,
                 VectorType& v,
                 T& b,
                 T& c,
                 const typename VectorType::size_type j = 0)
{
    v = u;

    FREEAML_ASSERT(j < v.size());

    const T norm = u.l2_norm();

    /* if u is the null vector, then it is already a multiple of e_j */
    if (norm == T{0})
    {
        b = T{0};
        c = T{1};
    }
    else
    {
        b = T{2};
        c = (u[j] >= T{0}) ? T{1} : T{-1};

        /* only the j-th element of v needs to be adjusted */
        v[j] += c * norm;
        v /= v.l2_norm();
    }
}

/**
 * @brief Computes the entries of the Givens rotation matrix for a 2d vector
 *        <tt>(x,y)</tt>.
 * @param x The @c x component of the 2d vector.
 * @param y The @c y component of the 2d vector.
 * @param c The first parameter which defines the Given rotation matrix for
 *        <tt>(x,y)</tt>.
 * @param s The second parameter which defines the Given rotation matrix for
 *        <tt>(x,y)</tt>.
 * @note The Givens rotation matrix has the form <tt>G = [[c -s, s c]]</tt>.
 */
template<typename T>
void givens(const T& x, const T& y, T& c, T& s)
{
    /* compute the magnitude of (x,y) */
    const T norm = std::sqrt(x * x + y * y);

    if (norm > T{0})
    {
        c = x / norm;
        s = -y / norm;
    }
    else
    {
        c = T{0};
        s = T{0};
    }
}

/**
 * @brief Applies a Givens rotation to a 2d vector <tt>(x,y)</tt>.
 * @param c The first parameter which defines the Given rotation matrix.
 * @param s The second parameter which defines the Given rotation matrix.
 * @param x The @c x component of the 2d vector.
 * @param y The @c y component of the 2d vector.
 * @note The Givens rotation matrix has the form <tt>G = [[c -s, s c]]</tt>.
 */
template<typename T>
void givens_rotation(const T& c, const T& s, T& x, T& y)
{
    const T _x = x;
    const T _y = y;

    x = c * _x - s * _y;
    y = s * _x + c * _y;
}

} /* namespace freeaml */
