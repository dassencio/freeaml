# target: run-tests executable based on Google Test
set(target run-tests)

# source files for generating the target
add_executable(
    ${target}
    BiconjugateGradientStabilized.cpp
    BidiagonalFactorization.cpp
    ConjugateGradient.cpp
    GaussianElimination.cpp
    GaussSeidel.cpp
    GeneralizedMinimumResidual.cpp
    IncompleteCholeskyConjugateGradient.cpp
    Jacobi.cpp
    Matrix.cpp
    MinimumResidual.cpp
    NormalEquations.cpp
    SparseMatrix.cpp
    SteepestDescent.cpp
    SuccessiveOverRelaxation.cpp
    Vector.cpp
    WeightedJacobi.cpp
)

# link the target to Google Test
target_link_libraries(${target} gtest_main)

# set the location of the header files used to build the target
target_include_directories(
    ${target}
    PUBLIC "${CMAKE_SOURCE_DIR}/include"                     # freeAML
    PUBLIC "${CMAKE_SOURCE_DIR}/contrib/googletest/include"  # Google Test
)

# register the target with ctest
add_test(NAME ${target} COMMAND ${target})

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
    "${CMAKE_CXX_FLAGS} -DFREEAML_DEBUG -O0 -Wall -Wextra -Werror"
)
