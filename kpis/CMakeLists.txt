# target: run-kpis executable based on Google Test
set(target run-kpis)

# source files for generating the target
add_executable(
    ${target}
    Matrix.cpp
    Vector.cpp
)

# link the target to Google Test
target_link_libraries(${target} gtest_main)

# set the location of the header files used to build the target
target_include_directories(
    ${target}
    PUBLIC "${CMAKE_SOURCE_DIR}/kpis"                        # KPI headers
    PUBLIC "${CMAKE_SOURCE_DIR}/include"                     # freeAML
    PUBLIC "${CMAKE_SOURCE_DIR}/contrib/googletest/include"  # Google Test
)

# place the target on the build directory root
set_target_properties(
    ${target}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}"
)

# set the compilation flags for building the target
set_target_properties(
    ${target}
    PROPERTIES
    COMPILE_FLAGS
    "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -Werror"
)
