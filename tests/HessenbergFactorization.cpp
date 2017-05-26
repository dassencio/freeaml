#include <HessenbergFactorization.h>
#include <Matrix.h>
#include <SparseMatrix.h>
#include <Vector.h>
#include <gtest/gtest.h>

TEST(HessenbergFactorizationTest, Factorize4x4Matrix)
{
    /* matrix to be factorized */
    freeaml::SparseMatrix<double> A = {{0.00, 4.20, 4.25, 0.00},
                                       {2.73, 1.23, 3.89, 2.11},
                                       {4.63, 0.00, 0.00, 0.00},
                                       {3.70, 0.00, 0.31, 3.84}};

    freeaml::Matrix<double> Q, H;

    /* compute the Hessenberg factorization A = QHQ^t (always succeeds) */
    freeaml::HessenbergFactorization::factorize(A, Q, H);

    /* compute Q^t*Q */
    freeaml::Matrix<double> QtQ = Q.transpose() * Q;

    /* compute Q*H*Q^t */
    freeaml::Matrix<double> QHQt = Q * H * Q.transpose();

    /* I is the 4×4 identity matrix */
    freeaml::Matrix<double> I = freeaml::identity_matrix<double>(4);

    EXPECT_LE((A - QHQt).max_norm(), 1.e-10);
    EXPECT_LE((QtQ - I).max_norm(), 1.e-10);
}

TEST(HessenbergFactorizationTest, Factorize4x4MatrixWithoutQ)
{
    /* matrix to be factorized */
    freeaml::SparseMatrix<double> A = {{0.00, 4.20, 4.25, 0.00},
                                       {2.73, 1.23, 3.89, 2.11},
                                       {4.63, 0.00, 0.00, 0.00},
                                       {3.70, 0.00, 0.31, 3.84}};

    freeaml::Matrix<double> H;

    /* compute the Hessenberg matrix H for A (always succeeds) */
    freeaml::HessenbergFactorization::factorize(A, H);

    /* exact Hessenberg matrix for A */
    freeaml::Matrix<double> Y = {{0.00, 4.77, 0.85, -3.49},
                                 {6.52, 3.22, -0.20, 1.55},
                                 {0.00, 3.91, -0.73, 0.15},
                                 {0.00, 0.00, 0.73, 2.57}};

    EXPECT_LE((H - Y).max_norm(), 1.e-2);
}

TEST(HessenbergFactorizationTest, Factorize5x5Matrix)
{
    /* matrix to be factorized */
    freeaml::Matrix<double> A = {{3.11, 4.20, 4.25, 1.73, 1.94},
                                 {2.73, 1.23, 3.89, 2.11, 0.48},
                                 {4.63, 0.00, 1.51, 2.67, 2.21},
                                 {3.70, 0.00, 0.31, 3.84, 1.44},
                                 {2.35, 0.67, 1.48, 4.82, 2.37}};

    freeaml::Matrix<double> Q, H;

    /* compute the Hessenberg factorization A = QHQ^t (always succeeds) */
    freeaml::HessenbergFactorization::factorize(A, Q, H);

    /* compute Q^t*Q */
    freeaml::Matrix<double> QtQ = Q.transpose() * Q;

    /* compute Q*H*Q^t */
    freeaml::Matrix<double> QHQt = Q * H * Q.transpose();

    /* I is the 5×5 identity matrix */
    freeaml::Matrix<double> I = freeaml::identity_matrix<double>(5);

    EXPECT_LE((A - QHQt).max_norm(), 1.e-10);
    EXPECT_LE((QtQ - I).max_norm(), 1.e-10);
}
