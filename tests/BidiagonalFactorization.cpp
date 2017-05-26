#include <BidiagonalFactorization.h>
#include <Matrix.h>
#include <SparseMatrix.h>
#include <Vector.h>
#include <gtest/gtest.h>

TEST(BidiagonalFactorizationTest, Factorize4x4Matrix)
{
    /* matrix to be factorized */
    freeaml::SparseMatrix<double> A = {{3.11, 4.20, 4.25, 1.73},
                                       {2.73, 1.23, 3.89, 2.11},
                                       {4.63, 0.00, 1.51, 2.67},
                                       {3.70, 0.00, 0.31, 3.84}};

    freeaml::Matrix<double> U, V;
    freeaml::SparseMatrix<double> Z;

    /* compute the bidiagonal factorization A = UZV^t */
    bool status = freeaml::BidiagonalFactorization::factorize(A, U, Z, V);

    EXPECT_TRUE(status);

    /* compute UZV^t */
    freeaml::Matrix<double> UZVt = U * Z * V.transpose();

    /* compute U^t*U and V^t*V */
    freeaml::Matrix<double> UtU = U.transpose() * U;
    freeaml::Matrix<double> VtV = V.transpose() * V;

    /* I is the 4×4 identity matrix */
    freeaml::Matrix<double> I = freeaml::identity_matrix<double>(4);

    EXPECT_LE((A - UZVt).max_norm(), 1.e-10);
    EXPECT_LE((UtU - I).max_norm(), 1.e-10);
    EXPECT_LE((VtV - I).max_norm(), 1.e-10);
}

TEST(BidiagonalFactorizationTest, Factorize6x4Matrix)
{
    /* matrix to be factorized */
    freeaml::Matrix<double> A = {{3.11, 4.20, 4.25, 1.73},
                                 {2.73, 1.23, 3.89, 2.11},
                                 {4.63, 0.00, 1.51, 2.67},
                                 {3.70, 0.00, 0.31, 3.84},
                                 {2.35, 0.67, 1.48, 4.82},
                                 {1.94, 0.48, 0.00, 1.44}};

    freeaml::Matrix<double> U, V;
    freeaml::SparseMatrix<double> Z;

    /* compute the bidiagonal factorization A = UZV^t */
    bool status = freeaml::BidiagonalFactorization::factorize(A, U, Z, V);

    EXPECT_TRUE(status);

    /* compute UZV^t */
    freeaml::Matrix<double> UZVt = U * Z * V.transpose();

    /* compute U^t*U and V^t*V */
    freeaml::Matrix<double> UtU = U.transpose() * U;
    freeaml::Matrix<double> VtV = V.transpose() * V;

    /* I4 and I6 are the 4×4 and 6×6 identity matrices respectively */
    freeaml::Matrix<double> I4 = freeaml::identity_matrix<double>(4);
    freeaml::Matrix<double> I6 = freeaml::identity_matrix<double>(6);

    EXPECT_LE((A - UZVt).max_norm(), 1.e-10);
    EXPECT_LE((UtU - I6).max_norm(), 1.e-10);
    EXPECT_LE((VtV - I4).max_norm(), 1.e-10);
}

TEST(BidiagonalFactorizationTest, Factorize4x6Matrix)
{
    /* matrix to be factorized */
    freeaml::Matrix<double> A = {{2.02, 0.07, 3.56, 0.96, 3.09, 2.93},
                                 {0.23, 1.13, 1.04, 4.56, 0.77, 2.89},
                                 {0.51, 1.93, 2.06, 4.80, 0.68, 2.28},
                                 {0.28, 0.83, 4.07, 2.25, 3.83, 0.13}};

    freeaml::Matrix<double> U, V;
    freeaml::SparseMatrix<double> Z;

    /* compute the bidiagonal factorization A = UZV^t */
    bool status = freeaml::BidiagonalFactorization::factorize(A, U, Z, V);

    EXPECT_TRUE(status);

    /* compute UZV^t */
    freeaml::Matrix<double> UZVt = U * Z * V.transpose();

    /* compute U^t*U and V^t*V */
    freeaml::Matrix<double> UtU = U.transpose() * U;
    freeaml::Matrix<double> VtV = V.transpose() * V;

    /* I4 and I6 are the 4×4 and 6×6 identity matrices respectively */
    freeaml::Matrix<double> I4 = freeaml::identity_matrix<double>(4);
    freeaml::Matrix<double> I6 = freeaml::identity_matrix<double>(6);

    EXPECT_LE((A - UZVt).max_norm(), 1.e-10);
    EXPECT_LE((UtU - I4).max_norm(), 1.e-10);
    EXPECT_LE((VtV - I6).max_norm(), 1.e-10);
}
