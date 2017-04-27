#include <Matrix.h>
#include <gtest/gtest.h>
#include <kpi.h>

TEST(MatrixKPI, MatrixMultiplication)
{
    KPI_BEGIN();

    freeaml::Matrix<double> M1 = freeaml::random_matrix<double>(2000, 2000);
    freeaml::Matrix<double> M2 = freeaml::random_matrix<double>(2000, 2000);

    KPI_MEASURE(M1 * M2);
}

TEST(MatrixKPI, MatrixVectorMultiplication)
{
    KPI_BEGIN();

    freeaml::Matrix<double> M = freeaml::random_matrix<double>(10000, 10000);
    freeaml::Vector<double> v = freeaml::random_vector<double>(10000);

    KPI_MEASURE(M * v);
}

TEST(MatrixKPI, VectorMatrixMultiplication)
{
    KPI_BEGIN();

    freeaml::Matrix<double> M = freeaml::random_matrix<double>(10000, 10000);
    freeaml::Vector<double> v = freeaml::random_vector<double>(10000);

    KPI_MEASURE(v * M);
}

TEST(MatrixKPI, MatrixTranspose)
{
    KPI_BEGIN();

    freeaml::Matrix<double> M = freeaml::random_matrix<double>(10000, 10000);

    KPI_MEASURE(M.transpose());
}
