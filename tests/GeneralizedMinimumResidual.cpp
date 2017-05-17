#include <GeneralizedMinimumResidual.h>
#include <Matrix.h>
#include <SparseMatrix.h>
#include <Vector.h>
#include <gtest/gtest.h>

TEST(GeneralizedMinimumResidualTest, Solve4x4LinearSystem)
{
    /* linear system matrix */
    freeaml::SparseMatrix<double> A = {{3.11, 4.20, 4.25, 0.00},
                                       {0.00, 0.00, 0.00, 0.30},
                                       {4.63, 0.00, 1.51, 4.38},
                                       {3.70, 0.00, 0.00, 3.84}};

    /* solution vector (initially set to zero) */
    freeaml::Vector<double> x = {0.0, 0.0, 0.0, 0.0};

    /* right-hand side */
    freeaml::Vector<double> b = {7.26, 0.3, 7.5, 7.54};

    /*
     * linear system solver: generalized minimum residual
     * maximum number of iterations: 4
     * residual tolerance: 1e-5
     */
    freeaml::GeneralizedMinimumResidual<double> lss(4, 1e-5);

    /* solve the linear system */
    bool status = lss.solve(A, x, b);

    EXPECT_TRUE(status);

    double residual = (A * x - b).l2_norm();

    EXPECT_LE(residual, lss.residual_tolerance());
    EXPECT_LE(lss.num_iterations(), lss.max_iterations());

    /* exact solution */
    freeaml::Vector<double> y = {1.0, 2.0, -1.0, 1.0};

    EXPECT_LE((x - y).l2_norm(), 1.e-5);
}

TEST(GeneralizedMinimumResidualTest, Solve5x5LinearSystem)
{
    /* linear system matrix */
    freeaml::Matrix<double> A = {{2.78, 1.57, 2.68, 4.38, 3.65},
                                 {2.56, 1.94, 1.98, 1.75, 0.83},
                                 {1.80, 1.25, 2.38, 0.62, 1.31},
                                 {2.59, 2.13, 0.65, 1.67, 0.37},
                                 {1.48, 2.52, 2.34, 0.53, 4.59}};

    /* solution vector (initially set to zero) */
    freeaml::Vector<double> x = {0.0, 0.0, 0.0, 0.0, 0.0};

    /* right-hand side */
    freeaml::Vector<double> b = {30.98, 13.03, 9.02, 12.69, 24.13};

    /*
     * linear system solver: generalized minimum residual
     * maximum number of iterations: 5
     * residual tolerance: 1e-6
     */
    freeaml::GeneralizedMinimumResidual<double> lss(5, 1e-6);

    /* solve the linear system */
    bool status = lss.solve(A, x, b);

    EXPECT_TRUE(status);

    double residual = (A * x - b).l2_norm();

    EXPECT_LE(residual, lss.residual_tolerance());
    EXPECT_LE(lss.num_iterations(), lss.max_iterations());

    /* exact solution */
    freeaml::Vector<double> y = {1.0, 2.0, -1.0, 3.0, 4.0};

    EXPECT_LE((x - y).l2_norm(), 1.e-5);
}
