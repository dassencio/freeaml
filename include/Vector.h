#pragma once

#include <debug.h>
#include <random>
#include <vector>

namespace freeaml
{
/**
 * @brief @c Vector<T> is an extension of @c std::vector<T> for mathematical
 * applications.
 *
 * This class stores a sequence of elements of type @c T. It overloads the
 * addition (+), the subtraction (-) and the multiplication (*) operators for
 * supporting common vector operations such as vector addition, vector
 * multiplication by scalar and vector dot product.
 *
 * Functions for computing @a L<SUP>p</SUP> norms as well as other commonly-used
 * mathematical operations are also provided in the class.
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

    /** @brief Constructs a vector with no elements. */
    Vector() = default;

    /**
     * @brief Constructs a vector with the contents of an initializer list.
     * @param init An initializer list holding elements of type @c T.
     */
    Vector(std::initializer_list<T> init);

    /**
     * @brief Constructs a vector with default-initialized elements.
     * @param n The length of the vector.
     */
    Vector(size_type n);

    /**
     * @brief Constructs a vector with all elements initialized with a value.
     * @param n The length of the vector.
     * @param x The initializing value for every element of the vector.
     */
    Vector(size_type n, const T& x);

    /**
     * @brief Constructs a vector with the contents of a range.
     * @param first An iterator pointing to the first range element.
     * @param last An iterator pointing to one-past-the-last range element.
     */
    template<class InputIterator>
    Vector(InputIterator first, InputIterator last);

    /**
     * @brief Copy constructor.
     * @param v The vector from which all elements will be copied.
     */
    Vector(const Vector& v) = default;

    /**
     * @brief Move constructor.
     * @param v The vector from which all elements will be moved.
     */
    Vector(Vector&& v) = default;

    /** @brief Destructor. */
    ~Vector() = default;

    /**
     * @brief Copy-assignment operator.
     * @param v The vector from which all elements will be copied.
     * @return A reference to @c *this.
     */
    Vector& operator=(const Vector& v) = default;

    /**
     * @brief Move-assignment operator.
     * @param v The vector from which all elements will be moved.
     * @return A reference to @c *this.
     */
    Vector& operator=(Vector&& v) = default;

    /**
     * @brief Multiplies all elements of the vector by a scalar.
     * @param c A scalar.
     * @return A reference to @c *this.
     */
    Vector& operator*=(const T& c);

    /**
     * @brief Divides all elements of the vector by a scalar.
     * @param c A scalar.
     * @return A reference to @c *this.
     */
    Vector& operator/=(const T& c);

    /**
     * @brief Performs element-wise addition-assignment with another vector.
     * @param v A vector.
     * @return A reference to @c *this.
     */
    Vector& operator+=(const Vector& v);

    /**
     * @brief Performs element-wise subtraction-assignment with another vector.
     * @param v A vector.
     * @return A reference to @c *this.
     */
    Vector& operator-=(const Vector& v);

    /**
     * @brief Computes the @a L<SUP>1</SUP>-norm of the vector.
     * @return The L<SUP>1</SUP>-norm of the vector.
     */
    T l1_norm() const;

    /**
     * @brief Computes the @a L<SUP>2</SUP>-norm of the vector.
     * @return The L<SUP>2</SUP>-norm of the vector.
     */
    T l2_norm() const;

    /**
     * @brief Computes the @a L<SUP>p</SUP>-norm of the vector.
     * @param p A scalar defining the norm to compute.
     * @return The @a L<SUP>p</SUP>-norm of the vector.
     */
    T lp_norm(const T& p) const;

    /**
     * @brief Computes the @a L<SUP>∞</SUP>-norm of the vector.
     * @return The @a L<SUP>∞</SUP>-norm of the vector.
     */
    T linf_norm() const;

    /**
     * @brief Computes the sum of all elements of the vector.
     * @return The sum of all elements of the vector.
     */
    T sum() const;

    /**
     * @brief Computes the arithmetic mean of the elements of the vector.
     * @return The arithmetic mean of the elements of the vector.
     */
    T mean() const;

}; /* class Vector<T> */

/**
 * @brief Computes the multiplication of a vector by a scalar on the right.
 * @param v A vector.
 * @param c A scalar.
 * @return A copy of @c v with all elements multiplied by @c c.
 */
template<typename T>
Vector<T> operator*(const Vector<T>& v, const T& c);

/**
 * @brief Computes the multiplication of a vector by a scalar on the left.
 * @param c A scalar.
 * @param v A vector.
 * @return A copy of @c v with all elements multiplied by @c c.
 */
template<typename T>
Vector<T> operator*(const T& c, const Vector<T>& v);

/**
 * @brief Computes the dot product of two equally-sized vectors.
 * @param v1 A vector.
 * @param v2 A vector.
 * @return The dot product of @c v1 and @c v2.
 */
template<typename T>
T operator*(const Vector<T>& v1, const Vector<T>& v2);

/**
 * @brief Computes the division of a vector by a scalar.
 * @param v A vector.
 * @param c A scalar.
 * @return A copy of @c v with all elements divided by @c c.
 */
template<typename T>
Vector<T> operator/(const Vector<T>& v, const T& c);

/**
 * @brief Computes the vectorial addition of two equally-sized vectors.
 * @param v1 A vector.
 * @param v2 A vector.
 * @return The element-wise sum of @c v1 and @c v2.
 */
template<typename T>
Vector<T> operator+(const Vector<T>& v1, const Vector<T>& v2);

/**
 * @brief Computes the vectorial difference of two equally-sized vectors.
 * @param v1 A vector.
 * @param v2 A vector.
 * @return The element-wise difference between @c v1 and @c v2.
 */
template<typename T>
Vector<T> operator-(const Vector<T>& v1, const Vector<T>& v2);

/**
 * @brief Computes the element-wise negation of a vector.
 * @param v A vector.
 * @return The element-wise negation of @c v.
 */
template<typename T>
Vector<T> operator-(const Vector<T>& v);

/**
 * @brief Prints the elements of a vector to an output stream.
 * @param stream An output stream to which the vector elements will be printed.
 * @param v A vector whose elements will be printed.
 * @return A reference to @c stream.
 */
template<typename T>
std::ostream& operator<<(std::ostream& stream, const Vector<T>& v);

/**
 * @brief Generates a random vector with elements within a specified range.
 * @param n The size of the vector to generate.
 * @param lower_bound The lower bound for the sample interval.
 * @param upper_bound The upper bound for the sample interval.
 * @return A vector with @c n elements sampled uniformly from
 *         @c [lower_bound, upper_bound].
 * @note This function was designed to work only with primitive integer and
 *       floating-point types (e.g. @c int, @c float, @c double etc.).
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
Vector<T>& Vector<T>::operator/=(const T& c)
{
    for (T& x : *this)
    {
        x /= c;
    }

    return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& v)
{
    FREEAML_ASSERT((*this).size() == v.size());

    for (size_type i = 0; i < (*this).size(); ++i)
    {
        (*this)[i] += v[i];
    }

    return *this;
}

template<typename T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& v)
{
    FREEAML_ASSERT((*this).size() == v.size());

    for (size_type i = 0; i < (*this).size(); ++i)
    {
        (*this)[i] -= v[i];
    }

    return *this;
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
Vector<T> operator*(const Vector<T>& v, const T& c)
{
    Vector<T> result = v;
    return (result *= c);
}

template<typename T>
Vector<T> operator*(const T& c, const Vector<T>& v)
{
    return (v * c);
}

template<typename T>
T operator*(const Vector<T>& v1, const Vector<T>& v2)
{
    FREEAML_ASSERT(v1.size() == v2.size());

    using size_type = typename Vector<T>::size_type;

    T dot_product{};

#ifdef _OPENMP
#pragma omp parallel
    {
        T local_dot_product{};

#pragma omp for nowait
        for (size_type i = 0; i < v1.size(); ++i)
        {
            local_dot_product += v1[i] * v2[i];
        }

#pragma omp critical
        {
            dot_product += local_dot_product;
        }
    }
#else
    /* serial implementation */
    for (size_type i = 0; i < v1.size(); ++i)
    {
        dot_product += v1[i] * v2[i];
    }
#endif /* #ifdef _OPENMP */

    return dot_product;
}

template<typename T>
Vector<T> operator/(const Vector<T>& v, const T& c)
{
    Vector<T> result = v;
    return (result /= c);
}

template<typename T>
Vector<T> operator+(const Vector<T>& v1, const Vector<T>& v2)
{
    FREEAML_ASSERT(v1.size() == v2.size());

    Vector<T> result = v1;
    return (result += v2);
}

template<typename T>
Vector<T> operator-(const Vector<T>& v1, const Vector<T>& v2)
{
    FREEAML_ASSERT(v1.size() == v2.size());

    Vector<T> result = v1;
    return (result -= v2);
}

template<typename T>
Vector<T> operator-(const Vector<T>& v)
{
    Vector<T> result;
    result.reserve(v.size());

    for (const T& x : v)
    {
        result.push_back(-x);
    }

    return result;
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
