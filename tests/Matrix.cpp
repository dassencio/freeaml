#include <Matrix.h>
#include <gtest/gtest.h>

TEST(MatrixTest, DefaultConstructor)
{
    const freeaml::Matrix<int> M;

    EXPECT_EQ(M.num_rows(), 0u);
    EXPECT_EQ(M.num_cols(), 0u);
}

TEST(MatrixTest, ConstructorWithInitializerList)
{
    const freeaml::Matrix<int> M1 = {};

    EXPECT_EQ(M1.num_rows(), 0u);
    EXPECT_EQ(M1.num_cols(), 0u);

    /*
     * M2 = [[ 0 1 2 ]
     *       [ 3 4 5 ]]
     */
    const freeaml::Matrix<int> M2 = {{0, 1, 2}, {3, 4, 5}};

    EXPECT_EQ(M2.num_rows(), 2u);
    EXPECT_EQ(M2.num_cols(), 3u);

    /* first row = {0, 1, 2} */
    EXPECT_EQ(M2(0, 0), 0);
    EXPECT_EQ(M2(0, 1), 1);
    EXPECT_EQ(M2(0, 2), 2);

    /* second row = {3, 4, 5} */
    EXPECT_EQ(M2(1, 0), 3);
    EXPECT_EQ(M2(1, 1), 4);
    EXPECT_EQ(M2(1, 2), 5);
}

TEST(MatrixTest, ConstructorWithVector)
{
    const freeaml::Matrix<int> M1(0, 0, freeaml::Vector<int>{});

    EXPECT_EQ(M1.num_rows(), 0u);
    EXPECT_EQ(M1.num_cols(), 0u);

    freeaml::Vector<int> elements = {0, 1, 2, 3, 4, 5};

    /*
     * M2 = M3 = [[ 0 1 2 ]
     *            [ 3 4 5 ]]
     */
    const freeaml::Matrix<int> M2(2, 3, elements);
    const freeaml::Matrix<int> M3(2, 3, std::move(elements));

    EXPECT_EQ(M2.num_rows(), 2u);
    EXPECT_EQ(M2.num_cols(), 3u);

    /* first row = {0, 1, 2} */
    EXPECT_EQ(M2(0, 0), 0);
    EXPECT_EQ(M2(0, 1), 1);
    EXPECT_EQ(M2(0, 2), 2);

    /* second row = {3, 4, 5} */
    EXPECT_EQ(M2(1, 0), 3);
    EXPECT_EQ(M2(1, 1), 4);
    EXPECT_EQ(M2(1, 2), 5);

    EXPECT_EQ(M2, M3);
}

TEST(MatrixTest, ConstructorWithMatrixDimensions)
{
    using size_type = typename freeaml::Matrix<int>::size_type;

    const freeaml::Matrix<int> M(5, 7);

    EXPECT_EQ(M.num_rows(), 5u);
    EXPECT_EQ(M.num_cols(), 7u);

    for (size_type i = 0; i < M.num_rows(); ++i)
    {
        for (size_type j = 0; j < M.num_cols(); ++j)
        {
            EXPECT_EQ(M(i, j), 0);
        }
    }
}

TEST(MatrixTest, ConstructorWithMatrixDimensionsAndDefaultValue)
{
    using size_type = typename freeaml::Matrix<int>::size_type;

    const freeaml::Matrix<int> M(5, 7, 2);

    EXPECT_EQ(M.num_rows(), 5u);
    EXPECT_EQ(M.num_cols(), 7u);

    for (size_type i = 0; i < M.num_rows(); ++i)
    {
        for (size_type j = 0; j < M.num_cols(); ++j)
        {
            EXPECT_EQ(M(i, j), 2);
        }
    }
}

TEST(MatrixTest, CopyConstructor)
{
    const freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Matrix<int> M2 = M1;

    EXPECT_EQ(M1, M2);
}

TEST(MatrixTest, MoveConstructor)
{
    freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Matrix<int> M2 = std::move(M1);

    EXPECT_EQ(M2, freeaml::Matrix<int>({{0, 1, 2}, {3, 4, 5}}));
}

TEST(MatrixTest, CopyAssignment)
{
    const freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    freeaml::Matrix<int> M2;

    M2 = M1;

    EXPECT_EQ(M1, M2);
}

TEST(MatrixTest, MoveAssignment)
{
    freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    freeaml::Matrix<int> M2;

    M2 = std::move(M1);

    EXPECT_EQ(M2, freeaml::Matrix<int>({{0, 1, 2}, {3, 4, 5}}));
}

TEST(MatrixTest, EqualityComparison)
{
    freeaml::Matrix<int> M1 = {{0, 1}, {2, 3}};
    freeaml::Matrix<int> M2 = M1;

    EXPECT_TRUE(M1 == M2);
    EXPECT_FALSE(M1 != M2);

    M2 = {{1, 2}, {3, 4}};

    EXPECT_TRUE(M1 != M2);
    EXPECT_FALSE(M1 == M2);

    M2 = {{0, 1}, {2, 3}, {4, 5}};

    EXPECT_TRUE(M1 != M2);
    EXPECT_FALSE(M1 == M2);
}

TEST(MatrixTest, MultiplicationByScalar)
{
    const freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Matrix<int> M2 = {{0, 2, 4}, {6, 8, 10}}; /* 2 * M1 */

    EXPECT_EQ(2 * M1, M2);
    EXPECT_EQ(M1 * 2, M2);
}

TEST(MatrixTest, DivisionByScalar)
{
    const freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Matrix<int> M2 = {{0, 0, 1}, {1, 2, 2}}; /* M1 / 2 */

    EXPECT_EQ(M1 / 2, M2);
}

TEST(MatrixTest, MatrixAddition)
{
    const freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Matrix<int> M2 = {{1, 2, 3}, {4, 5, 6}};
    const freeaml::Matrix<int> M3 = {{1, 3, 5}, {7, 9, 11}}; /* M1 + M2 */

    /* matrix addition is commutative */
    EXPECT_EQ(M1 + M2, M3);
    EXPECT_EQ(M2 + M1, M3);
}

TEST(MatrixTest, MatrixSubtraction)
{
    const freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Matrix<int> M2 = {{5, 4, 3}, {2, 1, 0}};
    const freeaml::Matrix<int> M3 = {{-5, -3, -1}, {1, 3, 5}}; /* M1 - M2 */

    /* matrix subtraction is anticommutative */
    EXPECT_EQ(M1 - M2, M3);
    EXPECT_EQ(M2 - M1, -M3);
}

TEST(MatrixTest, MatrixNegation)
{
    const freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Matrix<int> M2 = {{0, -1, -2}, {-3, -4, -5}}; /* -M1 */

    EXPECT_EQ(M1, -M2);
    EXPECT_EQ(M2, -M1);
    EXPECT_EQ(M1, -(-M1));
}

TEST(MatrixTest, MatrixTranspose)
{
    const freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Matrix<int> M2 = {{0, 3}, {1, 4}, {2, 5}}; /* M1^t */

    EXPECT_EQ(M1.transpose(), M2);
    EXPECT_EQ(M2.transpose(), M1);
    EXPECT_EQ(M1.transpose().transpose(), M1);
}

TEST(MatrixTest, IsSquare)
{
    const freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};

    EXPECT_FALSE(M1.is_square());

    const freeaml::Matrix<int> M2 = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};

    EXPECT_TRUE(M2.is_square());
}

TEST(MatrixTest, IsSymmetric)
{
    const freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};

    EXPECT_FALSE(M1.is_symmetric());

    const freeaml::Matrix<int> M2 = M1 + M1.transpose();

    EXPECT_TRUE(M2.is_symmetric());
}

TEST(MatrixTest, Resize)
{
    freeaml::Matrix<int> M;

    M.resize(2, 3);

    EXPECT_EQ(M, freeaml::Matrix<int>({{0, 0, 0}, {0, 0, 0}}));

    M.resize(2, 1);

    EXPECT_EQ(M, freeaml::Matrix<int>({{0}, {0}}));

    M.resize(2, 3, 1);

    EXPECT_EQ(M, freeaml::Matrix<int>({{0, 1, 1}, {0, 1, 1}}));

    M.resize(3, 2, 2);

    EXPECT_EQ(M, freeaml::Matrix<int>({{0, 1}, {0, 1}, {2, 2}}));

    M.resize(0, 0);

    EXPECT_EQ(M, freeaml::Matrix<int>{});
}

TEST(MatrixTest, Clear)
{
    freeaml::Matrix<int> M = {{0, 1, 2}, {3, 4, 5}};

    M.clear();

    EXPECT_EQ(M, freeaml::Matrix<int>{});
}

TEST(MatrixTest, Flatten)
{
    const freeaml::Matrix<int> M = {{0, 1, 2}, {3, 4, 5}};

    EXPECT_EQ(M.flatten(), freeaml::Vector<int>({0, 1, 2, 3, 4, 5}));
}

TEST(MatrixTest, MatrixMultiplication)
{
    const freeaml::Matrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Matrix<int> M2 = {{0, 1}, {2, 3}, {4, 5}};
    const freeaml::Matrix<int> M3 = {{10, 13}, {28, 40}}; /* M1 * M2 */

    EXPECT_EQ(M1 * M2, M3);
}

/* vector-matrix multiplication (the vector is interpreted as a row) */
TEST(MatrixTest, VectorMatrixMultiplication)
{
    const freeaml::Matrix<int> M = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Vector<int> v1 = {1, 2};
    const freeaml::Vector<int> v2 = {6, 9, 12}; /* v1 * M */

    EXPECT_EQ(v1 * M, v2);
}

/* matrix-vector multiplication (the vector is interpreted as a column) */
TEST(MatrixTest, MatrixVectorMultiplication)
{
    const freeaml::Matrix<int> M = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Vector<int> v1 = {1, 2, 3};
    const freeaml::Vector<int> v2 = {8, 26}; /* M * v1 */

    EXPECT_EQ(M * v1, v2);
}

TEST(MatrixTest, RandomIntegerMatrixWithinDefaultRange)
{
    const freeaml::Matrix<int> M = freeaml::random_matrix<int>(3, 4);

    EXPECT_EQ(M.num_rows(), 3u);
    EXPECT_EQ(M.num_cols(), 4u);
    EXPECT_EQ(M.flatten().size(), 12u);

    for (const int x : M.flatten())
    {
        EXPECT_GE(x, 0);
        EXPECT_LE(x, 1);
    }
}

TEST(MatrixTest, RandomIntegerMatrixWithinSpecifiedRange)
{
    const freeaml::Matrix<int> M = freeaml::random_matrix<int>(3, 4, -5, 5);

    EXPECT_EQ(M.num_rows(), 3u);
    EXPECT_EQ(M.num_cols(), 4u);
    EXPECT_EQ(M.flatten().size(), 12u);

    for (const int x : M.flatten())
    {
        EXPECT_GE(x, -5);
        EXPECT_LE(x, 5);
    }
}

TEST(MatrixTest, RandomFloatingPointMatrixWithinDefaultRange)
{
    const freeaml::Matrix<double> M = freeaml::random_matrix<double>(3, 4);

    EXPECT_EQ(M.num_rows(), 3u);
    EXPECT_EQ(M.num_cols(), 4u);
    EXPECT_EQ(M.flatten().size(), 12u);

    for (const double x : M.flatten())
    {
        EXPECT_GE(x, 0.0);
        EXPECT_LE(x, 1.0);
    }
}

TEST(MatrixTest, RandomFloatingPointMatrixWithinSpecifiedRange)
{
    const freeaml::Matrix<double> M =
        freeaml::random_matrix<double>(3, 4, -5.0, 5.0);

    EXPECT_EQ(M.num_rows(), 3u);
    EXPECT_EQ(M.num_cols(), 4u);
    EXPECT_EQ(M.flatten().size(), 12u);

    for (const double x : M.flatten())
    {
        EXPECT_GE(x, -5.0);
        EXPECT_LE(x, 5.0);
    }
}
