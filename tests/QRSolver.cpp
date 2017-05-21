#include <Matrix.h>
#include <QRSolver.h>
#include <SparseMatrix.h>
#include <Vector.h>
#include <gtest/gtest.h>

TEST(QRSolverTest, Solve5x3LinearLeastSquaresProblem)
{
    /* linear least-squares problem matrix */
    freeaml::SparseMatrix<double> A = {{3.47, 2.71, 0.00},
                                       {0.00, 2.15, 4.06},
                                       {0.00, 3.92, 0.19},
                                       {2.72, 0.00, 0.00},
                                       {1.08, 0.00, 0.87}};

    /* solution vector (initially set to zero) */
    freeaml::Vector<double> x = {0.0, 0.0, 0.0};

    /* right-hand side */
    freeaml::Vector<double> b = {1.41, 1.40, 4.06, 0.27, 4.20};

    /* solve the linear least-squares problem using QR factorization */
    bool status = freeaml::QRSolver::solve(A, x, b);

    EXPECT_TRUE(status);

    /* exact solution */
    freeaml::Vector<double> y = {0.159330, 0.713098, 0.185450};

    EXPECT_LE((x - y).l2_norm(), 1.e-6);
}

TEST(QRSolverTest, Solve6x4LinearLeastSquaresProblem)
{
    /* linear least-squares problem matrix */
    freeaml::Matrix<double> A = {{1.62, 0.43, 2.70, 4.45},
                                 {0.17, 2.65, 2.61, 2.28},
                                 {2.10, 3.24, 4.82, 4.95},
                                 {3.86, 0.27, 0.42, 1.28},
                                 {1.93, 2.47, 2.87, 2.15},
                                 {3.44, 1.29, 1.32, 4.45}};

    /* solution vector (initially set to zero) */
    freeaml::Vector<double> x = {0.0, 0.0, 0.0, 0.0};

    /* right-hand side */
    freeaml::Vector<double> b = {2.78, 3.39, 0.12, 1.01, 0.43, 3.39};

    /* solve the linear least-squares problem using QR factorization */
    bool status = freeaml::QRSolver::solve(A, x, b);

    EXPECT_TRUE(status);

    /* exact solution */
    freeaml::Vector<double> y = {-0.259086, 0.898976, -1.406186, 1.239579};

    EXPECT_LE((x - y).l2_norm(), 1.e-6);
}
