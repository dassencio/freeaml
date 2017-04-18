#pragma once

#include <debug.h>
#include <random>
#include <vector>

namespace freeaml
{
/**
 * Vector<T> is an extension of std::vector<T> for mathematical applications.
 *
 * This class overloads the + operator for adding vectors as well as the *
 * operator for supporting multiplication by scalar (both on the left and on
 * the right) or for computing the dot product between two vectors. Functions
 * for computing norms as well as other commonly-needed mathematical operations
 * are also provided.
 *
 * Support for OpenMP was added to the functions and operators which showed a
 * significant speedup when implemented using multiple threads.
 */
template<typename T>
class Vector : public std::vector<T>
{
public:
    using BaseVector = std::vector<T>;
    using size_type = typename BaseVector::size_type;

    /** @brief default constructor */
    Vector() = default;

    /** @brief constructs the vector with the contents of init */
    Vector(std::initializer_list<T> init);

    /** @brief constructs the vector with n default-initialized elements */
    Vector(size_type n);

    /** @brief constructs the vector with n elements initialized to x */
    Vector(size_type n, const T& x);

    /** @brief constructs the vector with the contents on range [first, last) */
    template<class InputIterator>
    Vector(InputIterator first, InputIterator last);

    /** @brief constructs the vector as a copy of other */
    Vector(const Vector& other) = default;

    /** @brief constructs the vector by moving the contents of other to it */
    Vector(Vector&& other) = default;

    /** @brief destructor */
    ~Vector() = default;

    /** @brief copies the elements of other to the vector */
    Vector& operator=(const Vector& other) = default;

    /** @brief moves the contents of other to the vector */
    Vector& operator=(Vector&& other) = default;

    /** @brief multiplies all elements of the vector by a scalar c */
    Vector& operator*=(const T& c);

    /** @brief computes the right-multiplication of the vector by a scalar c */
    Vector operator*(const T& c) const;

    /** @brief divides all elements of the vector by a scalar c */
    Vector& operator/=(const T& c);

    /** @brief computes the division of the vector by a scalar c */
    Vector operator/(const T& c) const;

    /** @brief adds other to the vector (vector addition) */
    Vector& operator+=(const Vector& other);

    /** @brief computes the result of adding other to the vector */
    Vector operator+(const Vector& other) const;

    /** @brief subtracts other from the vector (vector subtraction) */
    Vector& operator-=(const Vector& other);

    /** @brief computes the result of subtracting other from the vector */
    Vector operator-(const Vector& other) const;

    /** @brief computes the negation of the vector */
    Vector operator-() const;

    /** @brief computes the dot product of the vector with other */
    T operator*(const Vector& other) const;

    /** @brief computes the L^1-norm of the vector */
    T l1_norm() const;

    /** @brief computes the L^2-norm of the vector */
    T l2_norm() const;

    /** @brief computes the L^p-norm of the vector, for a given p value */
    T lp_norm(const T& p) const;

    /** @brief computes the L^inf-norm (or maximum norm) of the vector */
    T linf_norm() const;

    /** @brief computes the sum of all elements of the vector */
    T sum() const;

    /** @brief computes the arithmetic mean of the elements of the vector */
    T mean() const;

}; /* class Vector<T> */

/** @brief computes the left-multiplication of a vector v by a scalar c */
template<typename T>
Vector<T> operator*(const T& c, const Vector<T>& v);

/** @brief prints the elements of a vector v on an output stream */
template<typename T>
std::ostream& operator<<(std::ostream& stream, const Vector<T>& v);

/**
 * @brief generates a random vector with n elements within a given range
 * @param n the size of the vector to generate
 * @param lower_bound the lower bound for the sample interval
 * @param upper_bound the upper bound for the sample interval
 * @return the generated vector with n elements in [lower_bound, upper_bound]
 * @note this function was designed to work only with primitive integer and
 *       floating-point types (e.g. int, int32_t, float, double etc.)
 */
template<typename T>
Vector<T> random_vector(typename Vector<T>::size_type n,
                        const T& lower_bound = T{0},
                        const T& upper_bound = T{1});

/*******************************************************************************
 *
 *    FUNCTION DEFINITIONS
 *
 ******************************************************************************/

template<typename T>
Vector<T>::Vector(const std::initializer_list<T> init) : BaseVector(init)
{
    /* nothing needs to be done here */
}

template<typename T>
Vector<T>::Vector(const size_type n) : BaseVector(n)
{
    /* nothing needs to be done here */
}

template<typename T>
Vector<T>::Vector(const size_type n, const T& x) : BaseVector(n, x)
{
    /* nothing needs to be done here */
}

template<typename T>
template<typename InputIterator>
Vector<T>::Vector(const InputIterator first, const InputIterator last)
    : BaseVector(first, last)
{
    /* nothing needs to be done here */
}

template<typename T>
Vector<T>& Vector<T>::operator*=(const T& c)
{
    for (T& x : *this)
    {
        x *= c;
    }

    return *this;
}

template<typename T>
Vector<T> Vector<T>::operator*(const T& c) const
{
    Vector<T> v = *this;
    return (v *= c);
}

template<typename T>
Vector<T>& Vector<T>::operator/=(const T& c)
{
    for (T& x : *this)
    {
        x /= c;
    }

    return *this;
}

template<typename T>
Vector<T> Vector<T>::operator/(const T& c) const
{
    Vector<T> v = *this;
    return (v /= c);
}

template<typename T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& other)
{
    FREEAML_ASSERT((*this).size() == other.size());

    for (size_type i = 0; i < (*this).size(); ++i)
    {
        (*this)[i] += other[i];
    }

    return *this;
}

template<typename T>
Vector<T> Vector<T>::operator+(const Vector<T>& other) const
{
    FREEAML_ASSERT((*this).size() == other.size());

    Vector<T> v = *this;
    return (v += other);
}

template<typename T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& other)
{
    FREEAML_ASSERT((*this).size() == other.size());

    for (size_type i = 0; i < (*this).size(); ++i)
    {
        (*this)[i] -= other[i];
    }

    return *this;
}

template<typename T>
Vector<T> Vector<T>::operator-(const Vector<T>& other) const
{
    FREEAML_ASSERT((*this).size() == other.size());

    Vector<T> v = *this;
    return (v -= other);
}

template<typename T>
Vector<T> Vector<T>::operator-() const
{
    Vector<T> v;
    v.reserve((*this).size());

    for (const T& x : *this)
    {
        v.push_back(-x);
    }

    return v;
}

template<typename T>
T Vector<T>::operator*(const Vector<T>& other) const
{
    FREEAML_ASSERT((*this).size() == other.size());

    T dot_product{};

#ifdef _OPENMP
#pragma omp parallel
    {
        T local_dot_product{};

#pragma omp for nowait
        for (size_type i = 0; i < (*this).size(); ++i)
        {
            local_dot_product += (*this)[i] * other[i];
        }

#pragma omp critical
        {
            dot_product += local_dot_product;
        }
    }
#else
    /* serial implementation */
    for (size_type i = 0; i < (*this).size(); ++i)
    {
        dot_product += (*this)[i] * other[i];
    }
#endif /* #ifdef _OPENMP */

    return dot_product;
}

template<typename T>
T Vector<T>::l1_norm() const
{
    T norm{};

#ifdef _OPENMP
#pragma omp parallel
    {
        T local_norm{};

#pragma omp for nowait
        for (size_type i = 0; i < (*this).size(); ++i)
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
    for (const T& x : *this)
    {
        norm += std::abs(x);
    }
#endif /* #ifdef _OPENMP */

    return norm;
}

template<typename T>
T Vector<T>::l2_norm() const
{
    T norm{};

#ifdef _OPENMP
#pragma omp parallel
    {
        T local_norm = T{};

#pragma omp for nowait
        for (size_type i = 0; i < (*this).size(); ++i)
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
    for (const T& x : *this)
    {
        norm += std::abs(x) * std::abs(x);
    }
#endif /* #ifdef _OPENMP */

    return std::sqrt(norm);
}

template<typename T>
T Vector<T>::lp_norm(const T& p) const
{
    T norm{};

#ifdef _OPENMP
#pragma omp parallel
    {
        T local_norm{};

#pragma omp for nowait
        for (size_type i = 0; i < (*this).size(); ++i)
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
    for (const T& x : *this)
    {
        norm += std::pow(std::abs(x), p);
    }
#endif /* #ifdef _OPENMP */

    return std::pow(norm, T{1} / p);
}

template<typename T>
T Vector<T>::linf_norm() const
{
    T norm{};

#ifdef _OPENMP
#pragma omp parallel
    {
        T local_norm{};

#pragma omp for nowait
        for (size_type i = 0; i < (*this).size(); ++i)
        {
            local_norm = std::max(local_norm, std::abs((*this)[i]));
        }

#pragma omp critical
        {
            norm = std::max(local_norm, norm);
        }
    }
#else
    /* serial implementation */
    for (const T& x : *this)
    {
        norm = std::max(norm, std::abs(x));
    }
#endif /* #ifdef _OPENMP */

    return norm;
}

template<typename T>
T Vector<T>::sum() const
{
    T sum{};

#ifdef _OPENMP
#pragma omp parallel
    {
        T local_sum{};

#pragma omp for nowait
        for (size_type i = 0; i < (*this).size(); ++i)
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
    for (const T& x : *this)
    {
        sum += x;
    }
#endif /* #ifdef _OPENMP */

    return sum;
}

template<typename T>
T Vector<T>::mean() const
{
    if ((*this).empty() == false)
    {
        return sum() / static_cast<T>((*this).size());
    }
    else
    {
        return T{};
    }
}

template<typename T>
Vector<T> operator*(const T& c, const Vector<T>& v)
{
    return (v * c);
}

template<typename T>
std::ostream& operator<<(std::ostream& stream, const Vector<T>& v)
{
    stream << "[";

    for (typename Vector<T>::size_type i = 0; i < v.size(); ++i)
    {
        stream << v[i] << (i + 1 < v.size() ? ", " : "");
    }

    stream << "]\n";

    return stream;
}

template<typename T>
Vector<T> random_vector(const typename Vector<T>::size_type n,
                        const T& lower_bound /* = T{0} */,
                        const T& upper_bound /* = T{1} */)
{
    using DistributionType =
        typename std::conditional<std::is_integral<T>::value,
                                  std::uniform_int_distribution<T>,
                                  std::uniform_real_distribution<T>>::type;

    std::random_device device;
    std::mt19937_64 generator(device());
    DistributionType distribution(lower_bound, upper_bound);

    Vector<T> v;
    v.reserve(n);

    while (v.size() < n)
    {
        v.push_back(distribution(generator));
    }

    return v;
}

} /* namespace freeaml */
