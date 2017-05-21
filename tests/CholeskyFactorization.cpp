#include <CholeskyFactorization.h>
#include <Matrix.h>
#include <SparseMatrix.h>
#include <Vector.h>
#include <gtest/gtest.h>

TEST(CholeskyFactorizationTest, Factorize4x4Matrix)
{
    /* matrix to be factorized */
    freeaml::SparseMatrix<double> A = {{10.0, -1.0, 2.0, 0.0},
                                       {-1.0, 11.0, -1.0, 3.0},
                                       {2.0, -1.0, 10.0, -1.0},
                                       {0.0, 3.0, -1.0, 8.0}};

    freeaml::Matrix<double> L;

    /* compute the Cholesky factorization A = LL^t */
    bool status = freeaml::CholeskyFactorization::factorize(A, L);

    EXPECT_TRUE(status);

    /* compute LL^t */
    freeaml::Matrix<double> LLt = L * L.transpose();

    EXPECT_LE((A - LLt).max_norm(), 1.e-10);
}

TEST(CholeskyFactorizationTest, Factorize5x5Matrix)
{
    /* matrix to be factorized */
    freeaml::Matrix<double> A = {{1.49, 1.35, 1.25, 0.43, 1.49},
                                 {1.35, 1.50, 1.34, 0.53, 1.76},
                                 {1.25, 1.34, 1.79, 0.64, 2.01},
                                 {0.43, 0.53, 0.64, 0.28, 0.76},
                                 {1.49, 1.76, 2.01, 0.76, 2.50}};

    freeaml::Matrix<double> L;

    /* compute the Cholesky factorization A = LL^t */
    bool status = freeaml::CholeskyFactorization::factorize(A, L);

    EXPECT_TRUE(status);

    /* compute LL^t */
    freeaml::Matrix<double> LLt = L * L.transpose();

    EXPECT_LE((A - LLt).max_norm(), 1.e-10);
}
