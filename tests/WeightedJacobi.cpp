#include <Matrix.h>
#include <SparseMatrix.h>
#include <Vector.h>
#include <WeightedJacobi.h>
#include <gtest/gtest.h>

TEST(WeightedJacobiTest, Solve4x4LinearSystem)
{
    /* linear system matrix */
    freeaml::SparseMatrix<double> A = {{10.0, -1.0, 2.0, 0.0},
                                       {-1.0, 11.0, -1.0, 3.0},
                                       {2.0, -1.0, 10.0, -1.0},
                                       {0.0, 3.0, -1.0, 8.0}};

    /* solution vector (initially set to zero) */
    freeaml::Vector<double> x = {0.0, 0.0, 0.0, 0.0};

    /* right-hand side */
    freeaml::Vector<double> b = {6.0, 25.0, -11.0, 15.0};

    /*
     * linear system solver: weighted Jacobi
     * maximum number of iterations: 50
     * residual tolerance: 1e-7
     * iteration parameter (weight): 0.5
     */
    freeaml::WeightedJacobi<double> lss(50, 1e-7, 0.5);

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

TEST(WeightedJacobiTest, Iterate4x4LinearSystem)
{
    /* linear system matrix */
    freeaml::Matrix<double> A = {{10.0, -1.0, 2.0, 0.0},
                                 {-1.0, 11.0, -1.0, 3.0},
                                 {2.0, -1.0, 10.0, -1.0},
                                 {0.0, 3.0, -1.0, 8.0}};

    /* solution vector (initially set to zero) */
    freeaml::Vector<double> x = {0.0, 0.0, 0.0, 0.0};

    /* right-hand side */
    freeaml::Vector<double> b = {6.0, 25.0, -11.0, 15.0};

    /* perform 30 weighted-Jacobi iterations on the linear system */
    freeaml::WeightedJacobi<double>::iterate(A, x, b, 30, 0.66);

    /* exact solution */
    freeaml::Vector<double> y = {1.0, 2.0, -1.0, 1.0};

    EXPECT_LE((x - y).l2_norm(), 1.e-7);
}

TEST(WeightedJacobiTest, Solve5x5LinearSystem)
{
    /* linear system matrix */
    freeaml::Matrix<double> A = {{5.85, 1.05, 1.80, 0.78, 1.16},
                                 {1.05, 6.52, 0.57, 1.09, 0.34},
                                 {1.80, 0.57, 7.47, 1.01, 1.31},
                                 {0.78, 1.09, 1.01, 4.07, 0.15},
                                 {1.16, 0.34, 1.31, 0.15, 4.37}};

    /* solution vector (initially set to zero) */
    freeaml::Vector<double> x = {0.0, 0.0, 0.0, 0.0, 0.0};

    /* right-hand side */
    freeaml::Vector<double> b = {13.13, 18.15, 3.74, 14.76, 18.46};

    /*
     * linear system solver: weighted Jacobi
     * maximum number of iterations: 100
     * residual tolerance: 1e-8
     * iteration parameter (weight): 0.4
     */
    freeaml::WeightedJacobi<double> lss(100, 1e-8, 0.4);

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
