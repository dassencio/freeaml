[![Build Status](https://travis-ci.org/dassencio/freeaml.svg?branch=master)](https://travis-ci.org/dassencio/freeaml)

Description
===========

[FreeAML](https://github.com/dassencio/freeaml) is an easy-to-use, parallel
applied mathematics library written in C++. It boasts native integration with
[OpenMP](https://en.wikipedia.org/wiki/OpenMP)
and can therefore automatically use multiple cores to accelerate its implemented
algorithms (e.g. matrix multiplication, vectorial dot product etc.) on
shared-memory architectures.

FreeAML offers basic data structures tailored for mathematical applications
(vectors and matrices) as well as:

- linear system solvers
- linear least-squares solvers
- matrix factorization algorithms.

Integrating freeAML into your project can be done in minutes since it is a
header-only library. All you need to do is download the sources and include the
desired header files into your project.


License
=======

All code from this project is licensed under the FreeBSD license. See the
[`LICENSE`](https://github.com/dassencio/freeaml/tree/master/LICENSE) file for
more information.


Available tools
===============

The best way to learn how to use freeAML is by looking at its unit tests. Each
unit test is written in a very simple manner and is a self-contained example
of how a certain tool (e.g. a linear system solver) can be used. The following
list shows the available data structures and algorithms in the library, with
links pointing to their associated unit tests:

**Data structures**

- [Matrix](https://github.com/dassencio/freeaml/tree/master/tests/Matrix.cpp)
- [Sparse matrix](https://github.com/dassencio/freeaml/tree/master/tests/SparseMatrix.cpp)
- [Vector](https://github.com/dassencio/freeaml/tree/master/tests/Vector.cpp)

**Linear system solvers**

- [Biconjugate gradient stabilized (BiCGSTAB)](https://github.com/dassencio/freeaml/tree/master/tests/BiconjugateGradientStabilized.cpp)
- [Conjugate gradient (CG)](https://github.com/dassencio/freeaml/tree/master/tests/ConjugateGradient.cpp)
- [Conjugate gradient with incomplete-Cholesky preconditioning](https://github.com/dassencio/freeaml/tree/master/tests/IncompleteCholeskyConjugateGradient.cpp)
- [Gaussian elimination](https://github.com/dassencio/freeaml/tree/master/tests/GaussianElimination.cpp)
- [Gauss-Seidel](https://github.com/dassencio/freeaml/tree/master/tests/GaussSeidel.cpp)
- [Generalized minimum residual (GMRES)](https://github.com/dassencio/freeaml/tree/master/tests/GeneralizedMinimumResidual.cpp)
- [Jacobi](https://github.com/dassencio/freeaml/tree/master/tests/Jacobi.cpp)
- [Minimum residual (MINRES)](https://github.com/dassencio/freeaml/tree/master/tests/MinimumResidual.cpp)
- [Steepest descent](https://github.com/dassencio/freeaml/tree/master/tests/SteepestDescent.cpp)
- [Successive over-relaxation](https://github.com/dassencio/freeaml/tree/master/tests/SuccessiveOverRelaxation.cpp)
- [Weighted Jacobi](https://github.com/dassencio/freeaml/tree/master/tests/WeightedJacobi.cpp)

**Linear least-squares solvers**

- [Normal equations](https://github.com/dassencio/freeaml/tree/master/tests/NormalEquations.cpp)
- [QR solver](https://github.com/dassencio/freeaml/tree/master/tests/QRSolver.cpp)

**Matrix factorization algorithms**

- [Bidiagonal factorization](https://github.com/dassencio/freeaml/tree/master/tests/BidiagonalFactorization.cpp)
- [Cholesky factorization](https://github.com/dassencio/freeaml/tree/master/tests/CholeskyFactorization.cpp)
- [Hessenberg factorization](https://github.com/dassencio/freeaml/tree/master/tests/HessenbergFactorization.cpp)
- [PLU factorization](https://github.com/dassencio/freeaml/tree/master/tests/PLUFactorization.cpp)
- [QR factorization](https://github.com/dassencio/freeaml/tree/master/tests/QRFactorization.cpp)


Installation and usage instructions
===================================

FreeAML can be downloaded with the following command:

    git clone https://github.com/dassencio/freeaml.git

No installation is required: users can simply include the desired header files
from freeAML in their program sources and start focusing on what matters
to them.

If you wish to benefit from the native parallelization with OpenMP, make sure
that it is installed on your system (if you are using a popular Linux
distribution such as Ubuntu, it should already be there). Also, don't forget
to tell your compiler to link your application to the OpenMP library
(for instance, with `gcc`, you must add the `-fopenmp` option).

To compile and run the unit tests, go to the directory where you cloned the
freeAML sources, create a directory called `build` there, enter it and then run
the following commands in sequence:

    cmake -DUSE_OPENMP=ON ..
    make
    ./run-tests

If you do not wish to use OpenMP, simply replace `ON` with `OFF` on the command
above.


Doxygen documentation
=====================

FreeAML is thoroughly documented using the Doxygen format. To generate the
documentation as an HTML document, go to the directory where you cloned the
freeAML sources, create a directory called `build` there, enter it and then run
the following commands in sequence:

    cmake ..
    make doc

This will create a subdirectory called `html` containing an `index.html` file.
Open it on a browser to navigate through the documentation of the library
classes and functions.


Contributors & contact information
==================================

Diego Assencio / [diego@assencio.com](mailto:diego@assencio.com)
