#include <Jacobi.h>
#include <Matrix.h>
#include <Vector.h>
#include <gtest/gtest.h>

TEST(JacobiTest, Solve4x4LinearSystem)
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
     * Linear system solver: Jacobi
     * Maximum number of iterations: 50
     * Residual tolerance required: 1e-7
     */
    freeaml::Jacobi<double> lss(50, 1e-7);

    EXPECT_TRUE(lss.solve(A, x, b));

    /* check if the residual tolerance required was really achieved */
    EXPECT_LE((A * x - b).l2_norm(), lss.residual_tolerance());

    /* check if the number of iterations performed did not exceed the limit */
    EXPECT_LE(lss.num_iterations(), lss.max_iterations());

    std::cout << "Solution obtained: " << x << "\n"
              << "Residual achieved: " << (A * x - b).l2_norm() << "\n"
              << "Iterations performed: " << lss.num_iterations() << "\n";
}

TEST(JacobiTest, Iterate4x4LinearSystem)
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

    /* perform 20 Jacobi iterations on the linear system */
    freeaml::Jacobi<double>::iterate(A, x, b, 20);

    std::cout << "Solution obtained: " << x << "\n"
              << "Residual achieved: " << (A * x - b).l2_norm() << "\n";
}

TEST(JacobiTest, Solve5x5LinearSystem)
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
     * Linear system solver: Jacobi
     * Maximum number of iterations: 100
     * Residual tolerance required: 1e-8
     */
    freeaml::Jacobi<double> lss(100, 1e-8);

    EXPECT_TRUE(lss.solve(A, x, b));

    /* check if the residual tolerance required was really achieved */
    EXPECT_LE((A * x - b).l2_norm(), lss.residual_tolerance());

    /* check if the number of iterations performed did not exceed the limit */
    EXPECT_LE(lss.num_iterations(), lss.max_iterations());

    std::cout << "Solution obtained: " << x << "\n"
              << "Residual achieved: " << (A * x - b).l2_norm() << "\n"
              << "Iterations performed: " << lss.num_iterations() << "\n";
}
