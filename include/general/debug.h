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
