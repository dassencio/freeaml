cmake_minimum_required(VERSION 3.1)

project(freeAML CXX)

# set the project version number
set(VERSION_MAJOR 2)
set(VERSION_MINOR 0)

# require C++11
set(CMAKE_CXX_STANDARD 11)

# by default, enable OpenMP when building
option(USE_OPENMP "Use OpenMP" ON)

# if OpenMP is enabled, add the necessary compilation flags
if(USE_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    else()
        message(FATAL_ERROR "OpenMP could not be found.")
    endif()
endif()

# check if Doxygen is installed
find_package(Doxygen)

# if Doxygen is installed
if(DOXYGEN_FOUND)

    # parse and copy the Doxyfile to the build directory
    set(doxyfile_src ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(doxyfile_dst ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${doxyfile_src} ${doxyfile_dst} @ONLY)

    # create a "doc" target for building the Doxygen-based documentation
    add_custom_target(
        doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile_dst}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen..."
        VERBATIM
    )
endif()

# enable testing (for this directory and below) with ctest
enable_testing()

add_subdirectory(contrib)
add_subdirectory(tests)
add_subdirectory(kpis)
