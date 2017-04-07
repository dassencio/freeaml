#ifndef _freeAML_debug_h_
#define _freeAML_debug_h_


#include <string>
#include <iostream>
#include <cstdlib>


/** @brief stops the program execution until the user presses <enter> */
#define FREEAML_PAUSE std::cout << "\n\npress <Enter> to continue..."; getchar();


/** @brief terminates the program */
#define FREEAML_ABORT(errormsg) \
	aml::abort(errormsg, __FUNCTION__, __FILE__, __LINE__)


/** @brief terminates the program if some condition fails */
#ifdef FREEAML_DEBUG
#define FREEAML_ASSERT(condition) \
	aml::assert(condition, #condition, __PRETTY_FUNCTION__, __FILE__, __LINE__)
#else
#define FREEAML_ASSERT(condition)
#endif


namespace aml
{


void abort (const std::string errormsg,
            const std::string function,
            const std::string file,
            const int line)
{
	std::cerr << "\n\n" << std::string(37,'=')
	          << " ERROR " << std::string(36,'=')	<< "\n\n"
	          << "\n  error     : " << errormsg	<< "\n"
	          << "\n  function  : " << function	<< "\n"
	          << "\n  file      : " << file		<< "\n"
	          << "\n  line	    : " << line		<< "\n\n\n"
	          << std::endl;

	exit(EXIT_FAILURE);
}


void assert(const bool result,
            const std::string condition,
            const std::string function,
            const std::string file,
            const int line)
{
	if (!result)
	{
		abort("condition failed: " + condition, function, file, line);
	}
}


} /* end of namespace aml */

#endif /* _freeAML_debug_h */
