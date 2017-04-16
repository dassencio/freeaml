cmake_minimum_required(VERSION 3.1)

project(tests CXX)

# require C++11
set(CMAKE_CXX_STANDARD 11)

# set the compilation flags for the build
set(
    CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -DFREEAML_DEBUG -O3 -Wall -Wextra -Werror"
)

# add the OpenMP compilation flags if compilation with OpenMP is requested
OPTION (USE_OPENMP "Use OpenMP" ON)
IF(USE_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()

# build the gtest and gtest_main libraries
add_subdirectory(googletest)

# enable testing (for this directory and below) with ctest
enable_testing()

# source files for generating the test executable (run-tests)
add_executable(
    run-tests
    Vector.cpp
)

# link the test executable to gtest_main
target_link_libraries(run-tests gtest_main)

# indicate the location of the necessary header files (freeAML + Google Test)
target_include_directories(run-tests PUBLIC ../include)
target_include_directories(run-tests PUBLIC googletest/include)

# register run-tests as a test binary for ctest
add_test(NAME run-tests COMMAND run-tests)
