#include <IncompleteCholeskyConjugateGradient.h>
#include <SparseMatrix.h>
#include <Vector.h>
#include <gtest/gtest.h>

TEST(IncompleteCholeskyConjugateGradientTest, Solve4x4LinearSystem)
{
    /* linear system matrix */
    freeaml::SparseMatrix<double> A = {{10.0, -1.0, 0.0, 0.0},
                                       {-1.0, 11.0, -1.0, 0.0},
                                       {0.0, -1.0, 10.0, -1.0},
                                       {0.0, 0.0, -1.0, 8.0}};

    /* solution vector (initially set to zero) */
    freeaml::Vector<double> x = {0.0, 0.0, 0.0, 0.0};

    /* right-hand side */
    freeaml::Vector<double> b = {8.0, 22.0, -13.0, 9.0};

    /*
     * linear system solver: conjugate gradient with incomplete-Cholesky
     *                       preconditioning
     * maximum number of iterations: 4
     * residual tolerance: 1e-7
     */
    freeaml::IncompleteCholeskyConjugateGradient<double> lss(4, 1e-7);

    /* solve the linear system */
    bool status = lss.solve(A, x, b);

    EXPECT_TRUE(status);

    double residual = (A * x - b).l2_norm();

    EXPECT_LE(residual, lss.residual_tolerance());
    EXPECT_LE(lss.num_iterations(), lss.max_iterations());

    /* exact solution */
    freeaml::Vector<double> y = {1.0, 2.0, -1.0, 1.0};

    EXPECT_LE((x - y).l2_norm(), 1.e-7);
}

TEST(IncompleteCholeskyConjugateGradientTest, Solve5x5LinearSystem)
{
    /* linear system matrix */
    freeaml::SparseMatrix<double> A = {{3.24, 1.52, 1.50, 0.00, 0.00},
                                       {1.52, 4.66, 2.00, 1.45, 0.00},
                                       {1.50, 2.00, 3.68, 1.46, 0.00},
                                       {0.00, 1.45, 1.46, 3.00, 1.81},
                                       {0.00, 0.00, 0.00, 1.81, 5.56}};

    /* solution vector (initially set to zero) */
    freeaml::Vector<double> x = {0.0, 0.0, 0.0, 0.0, 0.0};

    /* right-hand side */
    freeaml::Vector<double> b = {4.78, 13.19, 6.2, 17.68, 27.67};

    /*
     * linear system solver: conjugate gradient with incomplete-Cholesky
     *                       preconditioning
     * maximum number of iterations: 5
     * residual tolerance: 1e-8
     */
    freeaml::IncompleteCholeskyConjugateGradient<double> lss(5, 1e-8);

    /* solve the linear system */
    bool status = lss.solve(A, x, b);

    EXPECT_TRUE(status);

    double residual = (A * x - b).l2_norm();

    EXPECT_LE(residual, lss.residual_tolerance());
    EXPECT_LE(lss.num_iterations(), lss.max_iterations());

    /* exact solution */
    freeaml::Vector<double> y = {1.0, 2.0, -1.0, 3.0, 4.0};

    EXPECT_LE((x - y).l2_norm(), 1.e-7);
}
