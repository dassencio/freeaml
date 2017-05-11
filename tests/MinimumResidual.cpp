#include <Matrix.h>
#include <MinimumResidual.h>
#include <Vector.h>
#include <gtest/gtest.h>

TEST(MinimumResidualTest, Solve4x4LinearSystem)
{
    /* linear system matrix */
    freeaml::Matrix<double> A = {{4.41, 1.36, 1.43, 3.09},
                                 {1.36, 4.64, 2.36, 0.50},
                                 {1.43, 2.36, 3.23, 0.92},
                                 {3.09, 0.50, 0.92, 0.87}};

    /* solution vector (initially set to zero) */
    freeaml::Vector<double> x = {0.0, 0.0, 0.0, 0.0};

    /* right-hand side */
    freeaml::Vector<double> b = {8.79, 8.78, 3.84, 4.04};

    /*
     * linear system solver: minimum residual
     * maximum number of iterations: 4
     * residual tolerance: 1e-7
     */
    freeaml::MinimumResidual<double> lss(4, 1e-7);

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

TEST(MinimumResidualTest, Solve5x5LinearSystem)
{
    /* linear system matrix */
    freeaml::Matrix<double> A = {{5.04, 6.82, 9.67, 4.17, 6.40},
                                 {6.82, 2.39, 6.29, 2.44, 5.52},
                                 {9.67, 6.29, 5.42, 6.09, 7.73},
                                 {4.17, 2.44, 6.09, 5.99, 8.89},
                                 {6.40, 5.52, 7.73, 8.89, 3.51}};

    /* solution vector (initially set to zero) */
    freeaml::Vector<double> x = {0.0, 0.0, 0.0, 0.0, 0.0};

    /* right-hand side */
    freeaml::Vector<double> b = {47.12, 34.71, 66.02, 56.49, 50.42};

    /*
     * linear system solver: minimum residual
     * maximum number of iterations: 5
     * residual tolerance: 1e-8
     */
    freeaml::MinimumResidual<double> lss(5, 1e-8);

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
