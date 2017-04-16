#pragma once

#include <cstdlib>
#include <iostream>

/**
 * @brief terminates the program if a given condition fails
 * @param condition the condition which will be tested
 */
#ifdef FREEAML_DEBUG
#define FREEAML_ASSERT(condition) \
    if ((condition) == false) \
    { \
        std::cerr << "\n\n" \
                  << "---- ASSERTION FAILED ----\n\n" \
                  << "  condition : " << #condition << "\n" \
                  << "  function  : " << __func__ << "\n" \
                  << "  file      : " << __FILE__ << "\n" \
                  << "  line      : " << __LINE__ << "\n\n"; \
        exit(EXIT_FAILURE); \
    }
#else
#define FREEAML_ASSERT(condition)
#endif
