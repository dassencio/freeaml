#include <Matrix.h>
#include <QRFactorization.h>
#include <SparseMatrix.h>
#include <Vector.h>
#include <gtest/gtest.h>

TEST(QRFactorizationTest, Factorize4x4Matrix)
{
    /* matrix to be factorized */
    freeaml::SparseMatrix<double> A = {{3.11, 4.20, 0.00, 1.73},
                                       {2.73, 1.23, 3.89, 0.00},
                                       {4.63, 0.00, 1.51, 2.67},
                                       {3.70, 0.00, 0.31, 3.84}};

    freeaml::Matrix<double> Q, R;

    /* compute the QR factorization A = QR (always succeeds) */
    freeaml::QRFactorization::factorize(A, Q, R);

    /* compute Q^t*Q */
    freeaml::Matrix<double> QtQ = Q.transpose() * Q;

    /* I is the 4×4 identity matrix */
    freeaml::Matrix<double> I = freeaml::identity_matrix<double>(4);

    EXPECT_LE((A - Q * R).max_norm(), 1.e-10);
    EXPECT_LE((QtQ - I).max_norm(), 1.e-10);
}

TEST(QRFactorizationTest, Factorize6x4Matrix)
{
    /* matrix to be factorized */
    freeaml::Matrix<double> A = {{3.11, 4.20, 4.25, 1.73},
                                 {2.73, 1.23, 3.89, 2.11},
                                 {4.63, 0.00, 1.51, 2.67},
                                 {3.70, 0.00, 0.31, 3.84},
                                 {2.35, 0.67, 1.48, 4.82},
                                 {1.94, 0.48, 0.00, 1.44}};

    freeaml::Matrix<double> Q, R;

    /* compute the QR factorization A = QR (always succeeds) */
    freeaml::QRFactorization::factorize(A, Q, R);

    /* compute Q^t*Q */
    freeaml::Matrix<double> QtQ = Q.transpose() * Q;

    /* I is the 4×4 identity matrix */
    freeaml::Matrix<double> I = freeaml::identity_matrix<double>(4);

    EXPECT_LE((A - Q * R).max_norm(), 1.e-10);
    EXPECT_LE((QtQ - I).max_norm(), 1.e-10);
}
