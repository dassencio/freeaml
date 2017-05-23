#pragma once

#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

/**
 * @brief Macro which must be added at the beginning of every KPI test.
 * @note If OpenMP is enabled, this macro will make each KPI test run using
 *       every possible number of threads (defined by omp_get_max_threads()).
 */
#ifdef _OPENMP
#define KPI_BEGIN()                              \
    do                                           \
    {                                            \
        int max_threads = omp_get_max_threads(); \
        for (int i = 1; i <= max_threads; ++i)   \
        {                                        \
            omp_set_num_threads(i)
#else
#define KPI_BEGIN() \
    do              \
    {
#endif /* _OPENMP */

/**
 * @brief Computes and prints the execution time of an expression.
 * @param expression A C++ expression.
 */
#ifdef _OPENMP
#define KPI_MEASURE(expression)                                              \
    auto t1 = std::chrono::steady_clock::now();                              \
    const auto& result = (expression);                                       \
    const volatile void* unused = reinterpret_cast<const void*>(&result);    \
    (void)unused;                                                            \
    auto t2 = std::chrono::steady_clock::now();                              \
    auto dt = std::chrono::duration<double, std::nano>(t2 - t1).count();     \
    std::cout << "Time taken with " << i << " threads (ns): " << dt << "\n"; \
    }                                                                        \
    }                                                                        \
    while (0)
#else
#define KPI_MEASURE(expression)                                           \
    auto t1 = std::chrono::steady_clock::now();                           \
    const auto& result = (expression);                                    \
    const volatile void* unused = reinterpret_cast<const void*>(&result); \
    (void)unused;                                                         \
    auto t2 = std::chrono::steady_clock::now();                           \
    auto dt = std::chrono::duration<double, std::nano>(t2 - t1).count();  \
    std::cout << "Time taken (ns): " << dt << "\n";                       \
    }                                                                     \
    while (0)
#endif /* _OPENMP */
