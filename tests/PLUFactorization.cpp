#include <Matrix.h>
#include <PLUFactorization.h>
#include <SparseMatrix.h>
#include <Vector.h>
#include <gtest/gtest.h>

TEST(PLUFactorizationTest, Factorize4x4MatrixWithoutPartialPivoting)
{
    /* matrix to be factorized */
    freeaml::SparseMatrix<double> A = {{0.00, 4.20, 4.25, 1.73},
                                       {0.00, 1.23, 3.89, 2.11},
                                       {4.63, 0.00, 1.51, 2.67},
                                       {3.70, 0.00, 0.31, 3.84}};

    freeaml::SparseMatrix<double> P;
    freeaml::Matrix<double> L, U;

    /* PLU factorizer without partial pivoting */
    freeaml::PLUFactorization fac(false);

    /* compute the PLU factorization A = PLU */
    bool status = fac.factorize(A, P, L, U);

    EXPECT_TRUE(status);

    /* compute P^t*P */
    freeaml::SparseMatrix<double> PtP = P.transpose() * P;

    /* I is the 4×4 identity matrix */
    freeaml::SparseMatrix<double> I =
        freeaml::identity_sparse_matrix<double>(4);

    EXPECT_LE((A - P * L * U).max_norm(), 1.e-10);
    EXPECT_LE((PtP - I).max_norm(), 1.e-10);
    EXPECT_NE(P, I);
}

TEST(PLUFactorizationTest, Factorize5x5MatrixWithPartialPivoting)
{
    /* matrix to be factorized */
    freeaml::Matrix<double> A = {{0.00, 1.57, 2.68, 4.38, 3.65},
                                 {2.56, 1.94, 1.98, 1.75, 0.83},
                                 {1.80, 1.25, 2.38, 0.62, 1.31},
                                 {2.59, 2.13, 0.65, 1.67, 0.37},
                                 {1.48, 2.52, 2.34, 0.53, 4.59}};

    freeaml::SparseMatrix<double> P;
    freeaml::Matrix<double> L, U;

    /* PLU factorizer with partial pivoting */
    freeaml::PLUFactorization fac(true);

    /* compute the PLU factorization A = PLU */
    bool status = fac.factorize(A, P, L, U);

    EXPECT_TRUE(status);

    /* compute P^t*P */
    freeaml::SparseMatrix<double> PtP = P.transpose() * P;

    /* I is the 5×5 identity matrix */
    freeaml::SparseMatrix<double> I =
        freeaml::identity_sparse_matrix<double>(5);

    EXPECT_LE((A - P * L * U).max_norm(), 1.e-10);
    EXPECT_LE((PtP - I).max_norm(), 1.e-10);
    EXPECT_NE(P, I);
}
