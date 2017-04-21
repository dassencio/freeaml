# target: run-tests executable based on Google Test
set(target run-tests)

# source files for generating the test executable
add_executable(
    ${target}
    Vector.cpp
)

# link the test executable to Google Test
target_link_libraries(${target} gtest_main)

# set the location of the header files used to build the test executable
target_include_directories(
    ${target}
    PUBLIC "${CMAKE_SOURCE_DIR}/include"                     # freeAML
    PUBLIC "${CMAKE_SOURCE_DIR}/contrib/googletest/include"  # Google Test
)

# register the test executable with ctest
add_test(NAME ${target} COMMAND ${target})

# place the test executable on the build directory root
set_target_properties(
    ${target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}"
)

# set the compilation flags for building the test executable
set_target_properties(
    ${target}
    PROPERTIES
    COMPILE_FLAGS
    "${CMAKE_CXX_FLAGS} -DFREEAML_DEBUG -O0 -Wall -Wextra -Werror"
)
