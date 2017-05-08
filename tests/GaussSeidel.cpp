#include <GaussSeidel.h>
#include <Matrix.h>
#include <Vector.h>
#include <gtest/gtest.h>

TEST(GaussSeidelTest, Solve4x4LinearSystem)
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
     * Linear system solver: Gauss-Seidel
     * Maximum number of iterations: 15
     * Residual tolerance required: 1e-7
     */
    freeaml::GaussSeidel<double> lss(15, 1e-7);

    /* solve the linear system */
    bool status = lss.solve(A, x, b);

    EXPECT_TRUE(status);

    double residual = (A * x - b).l2_norm();

    EXPECT_LE(residual, lss.residual_tolerance());
    EXPECT_LE(lss.num_iterations(), lss.max_iterations());

    std::cout << "Solution obtained: " << x << "\n"
              << "Residual achieved: " << residual << "\n"
              << "Iterations performed: " << lss.num_iterations() << "\n";
}

TEST(GaussSeidelTest, Iterate4x4LinearSystem)
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

    /* perform 10 Gauss-Seidel iterations on the linear system */
    freeaml::GaussSeidel<double>::iterate(A, x, b, 10);

    std::cout << "Solution obtained: " << x << "\n"
              << "Residual achieved: " << (A * x - b).l2_norm() << "\n";
}

TEST(GaussSeidelTest, Solve5x5LinearSystem)
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
     * Linear system solver: Gauss-Seidel
     * Maximum number of iterations: 20
     * Residual tolerance required: 1e-8
     */
    freeaml::GaussSeidel<double> lss(20, 1e-8);

    /* solve the linear system */
    bool status = lss.solve(A, x, b);

    EXPECT_TRUE(status);

    double residual = (A * x - b).l2_norm();

    EXPECT_LE(residual, lss.residual_tolerance());
    EXPECT_LE(lss.num_iterations(), lss.max_iterations());

    std::cout << "Solution obtained: " << x << "\n"
              << "Residual achieved: " << residual << "\n"
              << "Iterations performed: " << lss.num_iterations() << "\n";
}
