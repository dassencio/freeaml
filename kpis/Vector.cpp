#include <Vector.h>
#include <gtest/gtest.h>
#include <kpi.h>

TEST(VectorKPI, DotProduct)
{
    KPI_BEGIN();

    const freeaml::Vector<double> v1 =
        freeaml::random_vector<double>(100000000);
    const freeaml::Vector<double> v2 =
        freeaml::random_vector<double>(100000000);

    KPI_MEASURE(v1 * v2);
}

TEST(VectorKPI, L1Norm)
{
    KPI_BEGIN();

    const freeaml::Vector<double> v = freeaml::random_vector<double>(100000000);

    KPI_MEASURE(v.l1_norm());
}

TEST(VectorKPI, L2Norm)
{
    KPI_BEGIN();

    const freeaml::Vector<double> v = freeaml::random_vector<double>(100000000);

    KPI_MEASURE(v.l2_norm());
}

TEST(VectorKPI, LpNorm)
{
    KPI_BEGIN();

    const freeaml::Vector<double> v = freeaml::random_vector<double>(10000000);

    KPI_MEASURE(v.lp_norm(3.14));
}

TEST(VectorKPI, LInfNorm)
{
    KPI_BEGIN();

    const freeaml::Vector<double> v = freeaml::random_vector<double>(100000000);

    KPI_MEASURE(v.linf_norm());
}

TEST(VectorKPI, Sum)
{
    KPI_BEGIN();

    const freeaml::Vector<double> v = freeaml::random_vector<double>(100000000);

    KPI_MEASURE(v.sum());
}
