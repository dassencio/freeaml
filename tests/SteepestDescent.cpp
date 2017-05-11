#include <SteepestDescent.h>
#include <Matrix.h>
#include <Vector.h>
#include <gtest/gtest.h>

TEST(SteepestDescentTest, Solve4x4LinearSystem)
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

    /*
     * linear system solver: steepest descent
     * maximum number of iterations: 25
     * residual tolerance: 1e-7
     */
    freeaml::SteepestDescent<double> lss(25, 1e-7);

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

TEST(SteepestDescentTest, Solve5x5LinearSystem)
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
     * linear system solver: steepest descent
     * maximum number of iterations: 40
     * residual tolerance: 1e-8
     */
    freeaml::SteepestDescent<double> lss(40, 1e-8);

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
