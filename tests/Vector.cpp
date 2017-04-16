#include <Vector.h>
#include <gtest/gtest.h>

TEST(VectorTest, DefaultConstructor)
{
    freeaml::Vector<int> v;

    ASSERT_TRUE(v.empty());
    ASSERT_EQ(v.size(), 0u);
}

TEST(VectorTest, ConstructorWithInitializerList)
{
    freeaml::Vector<int> v = {1, 2, 3, 4};

    ASSERT_EQ(v.size(), 4u);

    EXPECT_EQ(v[0], 1);
    EXPECT_EQ(v[1], 2);
    EXPECT_EQ(v[2], 3);
    EXPECT_EQ(v[3], 4);
}

TEST(VectorTest, ConstructorWithLength)
{
    freeaml::Vector<int> v(10);

    ASSERT_EQ(v.size(), 10u);
}

TEST(VectorTest, ConstructorWithLengthAndDefaultValue)
{
    freeaml::Vector<int> v(10, 5);

    ASSERT_EQ(v.size(), 10u);

    for (const int x : v)
    {
        ASSERT_EQ(x, 5);
    }
}

TEST(VectorTest, ConstructorFromIteratorRange)
{
    freeaml::Vector<int> v1 = {1, 2, 3, 4, 5};
    freeaml::Vector<int> v2(v1.begin(), v1.end());

    ASSERT_EQ(v1, v2);
}

TEST(VectorTest, CopyAndMoveConstructor)
{
    freeaml::Vector<int> v1 = {1, 2, 3, 4, 5};
    freeaml::Vector<int> v2 = v1;

    ASSERT_EQ(v1, v2);

    freeaml::Vector<int> v3 = std::move(v2);

    ASSERT_EQ(v1, v3);
}

TEST(VectorTest, CopyAndMoveAssignment)
{
    freeaml::Vector<int> v1 = {1, 2, 3, 4, 5};
    freeaml::Vector<int> v2;

    v2 = v1;

    ASSERT_EQ(v1, v2);

    freeaml::Vector<int> v3;

    v3 = std::move(v2);

    ASSERT_EQ(v1, v3);
}

TEST(VectorTest, MultiplicationByScalar)
{
    freeaml::Vector<int> v1 = {1, 2, 3, 4, 5};
    freeaml::Vector<int> v2 = {2, 4, 6, 8, 10}; /* 2 * v1 */

    /* multiplication by scalar on the left */
    ASSERT_EQ(2 * v1, v2);

    /* multiplication by scalar on the right */
    ASSERT_EQ(v1 * 2, v2);

    v1 *= 2;

    ASSERT_EQ(v1, v2);
}

TEST(VectorTest, DivisionByScalar)
{
    freeaml::Vector<int> v1 = {1, 2, 3, 4, 5};
    freeaml::Vector<int> v2 = {0, 1, 1, 2, 2}; /* v1 / 2 */

    ASSERT_EQ(v1 / 2, v2);

    v1 /= 2;

    ASSERT_EQ(v1, v2);
}

TEST(VectorTest, VectorAddition)
{
    freeaml::Vector<int> v1 = {1, 2, 3, 4, 5};
    freeaml::Vector<int> v2 = {2, 3, 4, 5, 6};
    freeaml::Vector<int> v3 = {3, 5, 7, 9, 11}; /* v1 + v2 */

    ASSERT_EQ(v1 + v2, v3);

    v1 += v2;

    ASSERT_EQ(v1, v3);
}

TEST(VectorTest, VectorSubtraction)
{
    freeaml::Vector<int> v1 = {1, 2, 3, 4, 5};
    freeaml::Vector<int> v2 = {5, 4, 3, 2, 1};
    freeaml::Vector<int> v3 = {-4, -2, 0, 2, 4}; /* v1 - v2 */

    ASSERT_EQ(v1 - v2, v3);

    v1 -= v2;

    ASSERT_EQ(v1, v3);
}

TEST(VectorTest, VectorNegation)
{
    freeaml::Vector<int> v1 = {1, 2, 3, 4, 5};
    freeaml::Vector<int> v2 = {-1, -2, -3, -4, -5}; /* -v1 */

    ASSERT_EQ(v1, -v2);
}

TEST(VectorTest, DotProduct)
{
    freeaml::Vector<int> v1 = {1, 2, 3, 4, 5};
    freeaml::Vector<int> v2 = {2, 3, 4, 5, 6};

    ASSERT_EQ(v1 * v2, std::inner_product(v1.begin(), v1.end(), v2.begin(), 0));
}

TEST(VectorTest, L1Norm)
{
    freeaml::Vector<int> v = {1, -2, 3, -4, 5};

    ASSERT_EQ(v.l1_norm(), 15);
}

TEST(VectorTest, L2Norm)
{
    freeaml::Vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0};

    ASSERT_LT(std::abs(v.l2_norm() - 7.4162), 0.0001);
}

TEST(VectorTest, LpNorm)
{
    freeaml::Vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0};

    ASSERT_LT(std::abs(v.lp_norm(2.5) - 6.5583), 0.0001);
}

TEST(VectorTest, LInfNorm)
{
    freeaml::Vector<double> v = {-1.0, 2.0, -3.0, 4.0, -5.0};

    ASSERT_EQ(v.linf_norm(), 5.0);
}

TEST(VectorTest, SumOfElements)
{
    freeaml::Vector<int> v = {1, 2, 3, 4, 5};

    ASSERT_EQ(v.sum(), 15);
}

TEST(VectorTest, MeanValue)
{
    freeaml::Vector<double> v = {1.0, 2.0, 3.0, 4.0};

    ASSERT_LT(std::abs(v.mean() - 2.5), 1.e-10);
}

TEST(VectorTest, CreateRandomIntegerVector)
{
    freeaml::Vector<int> v = freeaml::random_vector<int>(10);

    ASSERT_EQ(v.size(), 10u);

    for (const int x : v)
    {
        EXPECT_GE(x, 0);
        EXPECT_LE(x, 1);
    }
}

TEST(VectorTest, CreateRandomIntegerVectorWithinSpecifiedRange)
{
    freeaml::Vector<int> v = freeaml::random_vector<int>(10, -5, 5);

    ASSERT_EQ(v.size(), 10u);

    for (const int x : v)
    {
        EXPECT_GE(x, -5);
        EXPECT_LE(x, 5);
    }
}

TEST(VectorTest, CreateRandomFloatingPointVector)
{
    freeaml::Vector<double> v = freeaml::random_vector<double>(10);

    ASSERT_EQ(v.size(), 10u);

    for (const double x : v)
    {
        EXPECT_GE(x, 0.0);
        EXPECT_LE(x, 1.0);
    }
}

TEST(VectorTest, CreateRandomFloatingPointVectorWithinSpecifiedRange)
{
    freeaml::Vector<double> v = freeaml::random_vector<double>(10, -1.0, 1.0);

    ASSERT_EQ(v.size(), 10u);

    for (const double x : v)
    {
        EXPECT_GE(x, -1.0);
        EXPECT_LE(x, 1.0);
    }
}
