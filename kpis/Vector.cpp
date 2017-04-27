#include <Vector.h>
#include <gtest/gtest.h>
#include <kpi.h>

TEST(VectorKPI, DotProduct)
{
    KPI_BEGIN();

    freeaml::Vector<double> v = freeaml::random_vector<double>(100000000);

    KPI_MEASURE(v * v);
}

TEST(VectorKPI, L1Norm)
{
    KPI_BEGIN();

    freeaml::Vector<double> v = freeaml::random_vector<double>(100000000);

    KPI_MEASURE(v.l1_norm());
}

TEST(VectorKPI, L2Norm)
{
    KPI_BEGIN();

    freeaml::Vector<double> v = freeaml::random_vector<double>(100000000);

    KPI_MEASURE(v.l2_norm());
}

TEST(VectorKPI, LpNorm)
{
    KPI_BEGIN();

    freeaml::Vector<double> v = freeaml::random_vector<double>(10000000);

    KPI_MEASURE(v.lp_norm(3.14));
}

TEST(VectorKPI, LInfNorm)
{
    KPI_BEGIN();

    freeaml::Vector<double> v = freeaml::random_vector<double>(100000000);

    KPI_MEASURE(v.linf_norm());
}

TEST(VectorKPI, Sum)
{
    KPI_BEGIN();

    freeaml::Vector<double> v = freeaml::random_vector<double>(100000000);

    KPI_MEASURE(v.sum());
}
