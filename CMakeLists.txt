cmake_minimum_required(VERSION 3.1)

project(freeAML CXX)

# require C++11
set(CMAKE_CXX_STANDARD 11)

# by default, enable OpenMP when building
option(USE_OPENMP "Use OpenMP" ON)

# if OpenMP is enabled, add the necessary compilation flags
if(USE_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# enable testing (for this directory and below) with ctest
enable_testing()

add_subdirectory(contrib)
add_subdirectory(tests)
