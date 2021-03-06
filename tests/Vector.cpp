#include <Vector.h>
#include <gtest/gtest.h>

TEST(VectorTest, DefaultConstructor)
{
    const freeaml::Vector<int> v;

    EXPECT_EQ(v.size(), 0u);
}

TEST(VectorTest, ConstructorWithInitializerList)
{
    const freeaml::Vector<int> v = {0, 1, 2, 3, 4};

    EXPECT_EQ(v.size(), 5u);

    EXPECT_EQ(v[0], 0);
    EXPECT_EQ(v[1], 1);
    EXPECT_EQ(v[2], 2);
    EXPECT_EQ(v[3], 3);
    EXPECT_EQ(v[4], 4);
}

TEST(VectorTest, ConstructorWithLength)
{
    const freeaml::Vector<int> v(10);

    EXPECT_EQ(v.size(), 10u);

    for (const int x : v)
    {
        EXPECT_EQ(x, 0);
    }
}

TEST(VectorTest, ConstructorWithLengthAndDefaultValue)
{
    const freeaml::Vector<int> v(10, 5);

    EXPECT_EQ(v.size(), 10u);

    for (const int x : v)
    {
        EXPECT_EQ(x, 5);
    }
}

TEST(VectorTest, ConstructorWithIteratorRange)
{
    const freeaml::Vector<int> v1 = {0, 1, 2, 3, 4};
    const freeaml::Vector<int> v2(v1.begin(), v1.end());

    EXPECT_EQ(v1, v2);
}

TEST(VectorTest, CopyConstructor)
{
    const freeaml::Vector<int> v1 = {0, 1, 2, 3, 4};
    const freeaml::Vector<int> v2 = v1;

    EXPECT_EQ(v1, v2);
}

TEST(VectorTest, MoveConstructor)
{
    freeaml::Vector<int> v1 = {0, 1, 2, 3, 4};
    const freeaml::Vector<int> v2 = std::move(v1);

    EXPECT_EQ(v2, freeaml::Vector<int>({0, 1, 2, 3, 4}));
}

TEST(VectorTest, CopyAssignment)
{
    const freeaml::Vector<int> v1 = {0, 1, 2, 3, 4};
    freeaml::Vector<int> v2;

    v2 = v1;

    EXPECT_EQ(v1, v2);
}

TEST(VectorTest, MoveAssignment)
{
    freeaml::Vector<int> v1 = {0, 1, 2, 3, 4};
    freeaml::Vector<int> v2;

    v2 = std::move(v1);

    EXPECT_EQ(v2, freeaml::Vector<int>({0, 1, 2, 3, 4}));
}

TEST(VectorTest, MultiplicationByScalar)
{
    const freeaml::Vector<int> v1 = {0, 1, 2, 3, 4};
    const freeaml::Vector<int> v2 = {0, 2, 4, 6, 8}; /* 2 * v1 */

    EXPECT_EQ(2 * v1, v2);
    EXPECT_EQ(v1 * 2, v2);
}

TEST(VectorTest, DivisionByScalar)
{
    const freeaml::Vector<int> v1 = {0, 1, 2, 3, 4};
    const freeaml::Vector<int> v2 = {0, 0, 1, 1, 2}; /* v1 / 2 */

    EXPECT_EQ(v1 / 2, v2);
}

TEST(VectorTest, VectorAddition)
{
    const freeaml::Vector<int> v1 = {0, 1, 2, 3, 4};
    const freeaml::Vector<int> v2 = {1, 2, 3, 4, 5};
    const freeaml::Vector<int> v3 = {1, 3, 5, 7, 9}; /* v1 + v2 */

    /* vector addition is commutative */
    EXPECT_EQ(v1 + v2, v3);
    EXPECT_EQ(v2 + v1, v3);
}

TEST(VectorTest, VectorSubtraction)
{
    const freeaml::Vector<int> v1 = {0, 1, 2, 3, 4};
    const freeaml::Vector<int> v2 = {4, 3, 2, 1, 0};
    const freeaml::Vector<int> v3 = {-4, -2, 0, 2, 4}; /* v1 - v2 */

    /* vector subtraction is anticommutative */
    EXPECT_EQ(v1 - v2, v3);
    EXPECT_EQ(v2 - v1, -v3);
}

TEST(VectorTest, VectorNegation)
{
    const freeaml::Vector<int> v1 = {0, 1, 2, 3, 4};
    const freeaml::Vector<int> v2 = {0, -1, -2, -3, -4}; /* -v1 */

    EXPECT_EQ(v1, -v2);
    EXPECT_EQ(v2, -v1);
    EXPECT_EQ(v1, -(-v1));
}

TEST(VectorTest, DotProduct)
{
    const freeaml::Vector<int> v1 = {0, 1, 2, 3, 4};
    const freeaml::Vector<int> v2 = {4, 3, 2, 1, 0};

    /* dot product is commutative */
    EXPECT_EQ(v1 * v2, 10);
    EXPECT_EQ(v2 * v1, 10);
}

TEST(VectorTest, Norms)
{
    const freeaml::Vector<double> v1 = {0.0, 1.0, -2.0, 3.0, -4.0};
    const freeaml::Vector<double> v2 = -v1;

    EXPECT_LT(std::abs(v1.l1_norm() - 10.0000), 0.0001);
    EXPECT_LT(std::abs(v2.l1_norm() - 10.0000), 0.0001);

    EXPECT_LT(std::abs(v1.l2_norm() - 5.4772), 0.0001);
    EXPECT_LT(std::abs(v2.l2_norm() - 5.4772), 0.0001);

    EXPECT_LT(std::abs(v1.lp_norm(1.0) - 10.0000), 0.0001);
    EXPECT_LT(std::abs(v2.lp_norm(1.0) - 10.0000), 0.0001);

    EXPECT_LT(std::abs(v1.lp_norm(2.0) - 5.4772), 0.0001);
    EXPECT_LT(std::abs(v2.lp_norm(2.0) - 5.4772), 0.0001);

    EXPECT_LT(std::abs(v1.lp_norm(2.5) - 4.9402), 0.0001);
    EXPECT_LT(std::abs(v2.lp_norm(2.5) - 4.9402), 0.0001);

    EXPECT_LT(std::abs(v1.lp_norm(3.14) - 4.5816), 0.0001);
    EXPECT_LT(std::abs(v2.lp_norm(3.14) - 4.5816), 0.0001);

    EXPECT_EQ(v1.linf_norm(), 4.0);
    EXPECT_EQ(v2.linf_norm(), 4.0);
}

TEST(VectorTest, ElementSum)
{
    const freeaml::Vector<int> v = {0, 1, 2, 3, 4};

    EXPECT_EQ(v.sum(), 10);
}

TEST(VectorTest, MeanValue)
{
    const freeaml::Vector<double> v = {1.0, 2.0, 3.0, 4.0};

    ASSERT_LT(std::abs(v.mean() - 2.5), 1.e-10);
}

TEST(VectorTest, RandomIntegerVectorWithinDefaultRange)
{
    const freeaml::Vector<int> v = freeaml::random_vector<int>(10);

    EXPECT_EQ(v.size(), 10u);

    for (const int x : v)
    {
        EXPECT_GE(x, 0);
        EXPECT_LE(x, 1);
    }
}

TEST(VectorTest, RandomIntegerVectorWithinSpecifiedRange)
{
    const freeaml::Vector<int> v = freeaml::random_vector<int>(10, -5, 5);

    EXPECT_EQ(v.size(), 10u);

    for (const int x : v)
    {
        EXPECT_GE(x, -5);
        EXPECT_LE(x, 5);
    }
}

TEST(VectorTest, RandomFloatingPointVectorWithinDefaultRange)
{
    const freeaml::Vector<double> v = freeaml::random_vector<double>(10);

    EXPECT_EQ(v.size(), 10u);

    for (const double x : v)
    {
        EXPECT_GE(x, 0.0);
        EXPECT_LE(x, 1.0);
    }
}

TEST(VectorTest, RandomFloatingPointVectorWithinSpecifiedRange)
{
    const freeaml::Vector<double> v =
        freeaml::random_vector<double>(10, -5.0, 5.0);

    EXPECT_EQ(v.size(), 10u);

    for (const double x : v)
    {
        EXPECT_GE(x, -5.0);
        EXPECT_LE(x, 5.0);
    }
}
