#include <SparseMatrix.h>
#include <gtest/gtest.h>
#include <kpi.h>

TEST(SparseMatrixKPI, SparseMatrixSparseMatrixMultiplication)
{
    KPI_BEGIN();

    const freeaml::SparseMatrix<double> M1 =
        freeaml::random_sparse_matrix<double>(1000, 1000, 100000);
    const freeaml::SparseMatrix<double> M2 =
        freeaml::random_sparse_matrix<double>(1000, 1000, 100000);

    KPI_MEASURE(M1 * M2);
}

TEST(SparseMatrixKPI, SparseMatrixMatrixMultiplication)
{
    KPI_BEGIN();

    const freeaml::SparseMatrix<double> M1 =
        freeaml::random_sparse_matrix<double>(1000, 1000, 100000);
    const freeaml::Matrix<double> M2 =
        freeaml::random_matrix<double>(1000, 1000);

    KPI_MEASURE(M1 * M2);
}

TEST(SparseMatrixKPI, MatrixSparseMatrixMultiplication)
{
    KPI_BEGIN();

    const freeaml::Matrix<double> M1 =
        freeaml::random_matrix<double>(1000, 1000);
    const freeaml::SparseMatrix<double> M2 =
        freeaml::random_sparse_matrix<double>(1000, 1000, 100000);

    KPI_MEASURE(M1 * M2);
}

TEST(SparseMatrixKPI, SparseMatrixVectorMultiplication)
{
    KPI_BEGIN();

    const freeaml::SparseMatrix<double> M =
        freeaml::random_sparse_matrix<double>(30000, 30000, 10000000);
    const freeaml::Vector<double> v = freeaml::random_vector<double>(30000);

    KPI_MEASURE(M * v);
}

TEST(SparseMatrixKPI, VectorSparseMatrixMultiplication)
{
    KPI_BEGIN();

    const freeaml::SparseMatrix<double> M =
        freeaml::random_sparse_matrix<double>(30000, 30000, 10000000);
    const freeaml::Vector<double> v = freeaml::random_vector<double>(30000);

    KPI_MEASURE(v * M);
}

TEST(SparseMatrixKPI, SparseMatrixSparseMatrixAddition)
{
    KPI_BEGIN();

    const freeaml::SparseMatrix<double> M1 =
        freeaml::random_sparse_matrix<double>(30000, 30000, 10000000);
    const freeaml::SparseMatrix<double> M2 =
        freeaml::random_sparse_matrix<double>(30000, 30000, 10000000);

    KPI_MEASURE(M1 + M2);
}

TEST(SparseMatrixKPI, SparseMatrixSparseMatrixSubtraction)
{
    KPI_BEGIN();

    const freeaml::SparseMatrix<double> M1 =
        freeaml::random_sparse_matrix<double>(30000, 30000, 10000000);
    const freeaml::SparseMatrix<double> M2 =
        freeaml::random_sparse_matrix<double>(30000, 30000, 10000000);

    KPI_MEASURE(M1 - M2);
}
