#include <SparseMatrix.h>
#include <gtest/gtest.h>

TEST(SparseMatrixTest, DefaultConstructor)
{
    const freeaml::SparseMatrix<int> M;

    EXPECT_EQ(M.num_rows(), 0u);
    EXPECT_EQ(M.num_cols(), 0u);
    EXPECT_EQ(M.num_stored(), 0u);
}

TEST(SparseMatrixTest, ConstructorWithInitializerList)
{
    const freeaml::SparseMatrix<int> M1 = {};

    EXPECT_EQ(M1.num_rows(), 0u);
    EXPECT_EQ(M1.num_cols(), 0u);
    EXPECT_EQ(M1.num_stored(), 0u);

    /*
     * M2 = [[ 0 1 2 ]
     *       [ 3 4 5 ]]
     */
    const freeaml::SparseMatrix<int> M2 = {{0, 1, 2}, {3, 4, 5}};

    EXPECT_EQ(M2.num_rows(), 2u);
    EXPECT_EQ(M2.num_cols(), 3u);
    EXPECT_EQ(M2.num_stored(), 5u);

    /* first row = {0, 1, 2} */
    EXPECT_EQ(M2(0, 0), 0);
    EXPECT_EQ(M2(0, 1), 1);
    EXPECT_EQ(M2(0, 2), 2);

    /* second row = {3, 4, 5} */
    EXPECT_EQ(M2(1, 0), 3);
    EXPECT_EQ(M2(1, 1), 4);
    EXPECT_EQ(M2(1, 2), 5);
}

TEST(SparseMatrixTest, ConstructorWithMatrixDimensions)
{
    using size_type = typename freeaml::SparseMatrix<int>::size_type;

    const freeaml::SparseMatrix<int> M(5, 7);

    EXPECT_EQ(M.num_rows(), 5u);
    EXPECT_EQ(M.num_cols(), 7u);
    EXPECT_EQ(M.num_stored(), 0u);

    for (size_type i = 0; i < M.num_rows(); ++i)
    {
        for (size_type j = 0; j < M.num_cols(); ++j)
        {
            EXPECT_EQ(M(i, j), 0);
        }
    }
}

TEST(SparseMatrixTest, ConstructorWithVector)
{
    const freeaml::SparseMatrix<int> M1(0, 0, freeaml::Vector<int>{});

    EXPECT_EQ(M1.num_rows(), 0u);
    EXPECT_EQ(M1.num_cols(), 0u);
    EXPECT_EQ(M1.num_stored(), 0u);

    const freeaml::Vector<int> elements = {0, 1, 2, 3, 4, 5};

    /*
     * M2 = [[ 0 1 2 ]
     *       [ 3 4 5 ]]
     */
    const freeaml::SparseMatrix<int> M2(2, 3, elements);

    EXPECT_EQ(M2.num_rows(), 2u);
    EXPECT_EQ(M2.num_cols(), 3u);
    EXPECT_EQ(M2.num_stored(), 5u);

    /* first row = {0, 1, 2} */
    EXPECT_EQ(M2(0, 0), 0);
    EXPECT_EQ(M2(0, 1), 1);
    EXPECT_EQ(M2(0, 2), 2);

    /* second row = {3, 4, 5} */
    EXPECT_EQ(M2(1, 0), 3);
    EXPECT_EQ(M2(1, 1), 4);
    EXPECT_EQ(M2(1, 2), 5);
}

TEST(SparseMatrixTest, CopyConstructor)
{
    const freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::SparseMatrix<int> M2 = M1;

    EXPECT_EQ(M1.num_stored(), 5u);
    EXPECT_EQ(M2.num_stored(), 5u);
    EXPECT_EQ(M1, M2);
}

TEST(SparseMatrixTest, MoveConstructor)
{
    freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::SparseMatrix<int> M2 = std::move(M1);

    EXPECT_EQ(M2.num_stored(), 5u);
    EXPECT_EQ(M2, freeaml::SparseMatrix<int>({{0, 1, 2}, {3, 4, 5}}));
}

TEST(SparseMatrixTest, CopyAssignment)
{
    const freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    freeaml::SparseMatrix<int> M2;

    M2 = M1;

    EXPECT_EQ(M1.num_stored(), 5u);
    EXPECT_EQ(M2.num_stored(), 5u);
    EXPECT_EQ(M1, M2);
}

TEST(SparseMatrixTest, MoveAssignment)
{
    freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    freeaml::SparseMatrix<int> M2;

    M2 = std::move(M1);

    EXPECT_EQ(M2.num_stored(), 5u);
    EXPECT_EQ(M2, freeaml::SparseMatrix<int>({{0, 1, 2}, {3, 4, 5}}));
}

TEST(SparseMatrixTest, EqualityComparison)
{
    freeaml::SparseMatrix<int> M1 = {{0, 1}, {2, 3}};
    freeaml::SparseMatrix<int> M2 = M1;

    EXPECT_TRUE(M1 == M2);
    EXPECT_FALSE(M1 != M2);

    M2 = {{1, 2}, {3, 4}};

    EXPECT_TRUE(M1 != M2);
    EXPECT_FALSE(M1 == M2);

    M2 = {{0, 1}, {2, 3}, {4, 5}};

    EXPECT_TRUE(M1 != M2);
    EXPECT_FALSE(M1 == M2);
}

TEST(SparseMatrixTest, MultiplicationByScalar)
{
    const freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::SparseMatrix<int> M2 = {{0, 2, 4}, {6, 8, 10}}; /* 2 * M1 */

    EXPECT_EQ(2 * M1, M2);
    EXPECT_EQ(M1 * 2, M2);
}

TEST(SparseMatrixTest, DivisionByScalar)
{
    const freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::SparseMatrix<int> M2 = {{0, 0, 1}, {1, 2, 2}}; /* M1 / 2 */

    EXPECT_EQ(M1 / 2, M2);
}

TEST(SparseMatrixTest, MatrixAddition)
{
    const freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::SparseMatrix<int> M2 = {{1, 2, 3}, {4, 5, 6}};
    const freeaml::SparseMatrix<int> M3 = {{1, 3, 5}, {7, 9, 11}}; /* M1 + M2 */

    /* matrix addition is commutative */
    EXPECT_EQ(M1 + M2, M3);
    EXPECT_EQ(M2 + M1, M3);

    /* matrix addition also works with dense matrices */
    const freeaml::Matrix<int> M4 = {{1, 2, 3}, {4, 5, 6}};
    const freeaml::Matrix<int> M5 = {{1, 3, 5}, {7, 9, 11}}; /* M1 + M4 */

    /* matrix addition is commutative */
    EXPECT_EQ(M1 + M4, M5);
    EXPECT_EQ(M4 + M1, M5);
}

TEST(SparseMatrixTest, MatrixSubtraction)
{
    const freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::SparseMatrix<int> M2 = {{5, 4, 3}, {2, 1, 0}};
    const freeaml::SparseMatrix<int> M3 = {{-5, -3, -1},
                                           {1, 3, 5}}; /* M1 - M2 */

    /* matrix subtraction is anticommutative */
    EXPECT_EQ(M1 - M2, M3);
    EXPECT_EQ(M2 - M1, -M3);

    /* matrix subtraction also works with dense matrices */
    const freeaml::Matrix<int> M4 = {{5, 4, 3}, {2, 1, 0}};
    const freeaml::Matrix<int> M5 = {{-5, -3, -1}, {1, 3, 5}}; /* M1 - M4 */

    /* matrix subtraction is anticommutative */
    EXPECT_EQ(M1 - M4, M5);
    EXPECT_EQ(M4 - M1, -M5);
}

TEST(SparseMatrixTest, MatrixNegation)
{
    const freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::SparseMatrix<int> M2 = {{0, -1, -2}, {-3, -4, -5}}; /* -M1 */

    EXPECT_EQ(M1, -M2);
    EXPECT_EQ(M2, -M1);
    EXPECT_EQ(M1, -(-M1));
}

TEST(SparseMatrixTest, MatrixTranspose)
{
    const freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::SparseMatrix<int> M2 = {{0, 3}, {1, 4}, {2, 5}}; /* M1^t */

    EXPECT_EQ(M1.transpose(), M2);
    EXPECT_EQ(M2.transpose(), M1);
    EXPECT_EQ(M1.transpose().transpose(), M1);
}

TEST(SparseMatrixTest, MaxNorm)
{
    const freeaml::SparseMatrix<int> M = {{0, -1, 2}, {-5, 4, -3}};

    EXPECT_EQ(M.max_norm(), 5);
}

TEST(SparseMatrixTest, IsSquare)
{
    const freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};

    EXPECT_FALSE(M1.is_square());

    const freeaml::SparseMatrix<int> M2 = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};

    EXPECT_TRUE(M2.is_square());
}

TEST(SparseMatrixTest, IsSymmetric)
{
    const freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};

    EXPECT_FALSE(M1.is_symmetric());

    const freeaml::SparseMatrix<int> M2 = M1 + M1.transpose();

    EXPECT_TRUE(M2.is_symmetric());
}

TEST(SparseMatrixTest, Resize)
{
    freeaml::SparseMatrix<int> M;

    M.resize(2, 3);

    EXPECT_EQ(M.num_stored(), 0u);
    EXPECT_EQ(M, freeaml::SparseMatrix<int>({{0, 0, 0}, {0, 0, 0}}));

    M(0, 1) = 1;
    M(1, 2) = 1;

    EXPECT_EQ(M.num_stored(), 2u);
    EXPECT_EQ(M, freeaml::SparseMatrix<int>({{0, 1, 0}, {0, 0, 1}}));

    M.resize(3, 2);

    EXPECT_EQ(M.num_stored(), 1u);
    EXPECT_EQ(M, freeaml::SparseMatrix<int>({{0, 1}, {0, 0}, {0, 0}}));

    M.resize(0, 0);

    EXPECT_EQ(M.num_stored(), 0u);
    EXPECT_EQ(M, freeaml::SparseMatrix<int>{});
}

TEST(SparseMatrixTest, Clear)
{
    freeaml::SparseMatrix<int> M = {{0, 1, 2}, {3, 4, 5}};

    M.clear();

    EXPECT_EQ(M.num_stored(), 0u);
    EXPECT_EQ(M, freeaml::SparseMatrix<int>{});
}

TEST(SparseMatrixTest, Flatten)
{
    const freeaml::SparseMatrix<int> M = {{0, 1, 2}, {3, 4, 5}};

    EXPECT_EQ(M.flatten(), freeaml::Vector<int>({0, 1, 2, 3, 4, 5}));
}

TEST(SparseMatrixTest, NumStored)
{
    freeaml::SparseMatrix<int> M = {{0, 1, 0}, {1, 0, 1}};

    EXPECT_EQ(M.num_stored(), 3u);

    /*
     * accessing a zero element in a non-const matrix will cause it to be
     * explicitly stored in the matrix with value equal to zero (bad!)
     */
    M(0, 0);

    EXPECT_EQ(M.num_stored(), 4u);
}

TEST(SparseMatrixTest, MatrixMultiplication)
{
    const freeaml::SparseMatrix<int> M1 = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::SparseMatrix<int> M2 = {{0, 1}, {2, 3}, {4, 5}};
    const freeaml::SparseMatrix<int> M3 = {{10, 13}, {28, 40}}; /* M1 * M2 */

    EXPECT_EQ(M1 * M2, M3);

    /* matrix multiplication also works with dense matrices */
    const freeaml::Matrix<int> M4 = {{0, 1}, {2, 3}, {4, 5}};
    const freeaml::Matrix<int> M5 = {{10, 13}, {28, 40}}; /* M1 * M4 */
    const freeaml::Matrix<int> M6 = {
        {3, 4, 5}, {9, 14, 19}, {15, 24, 33}}; /* M4 * M1 */

    /* (sparse matrix) * (dense matrix) = (dense matrix) */
    EXPECT_EQ(M1 * M4, M5);

    /* (dense matrix) * (sparse matrix) = (dense matrix) */
    EXPECT_EQ(M4 * M1, M6);
}

/* vector-matrix multiplication (the vector is interpreted as a row) */
TEST(SparseMatrixTest, VectorMatrixMultiplication)
{
    const freeaml::SparseMatrix<int> M = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Vector<int> v1 = {1, 2};
    const freeaml::Vector<int> v2 = {6, 9, 12}; /* v1 * M */

    EXPECT_EQ(v1 * M, v2);
}

/* matrix-vector multiplication (the vector is interpreted as a column) */
TEST(SparseMatrixTest, MatrixVectorMultiplication)
{
    const freeaml::SparseMatrix<int> M = {{0, 1, 2}, {3, 4, 5}};
    const freeaml::Vector<int> v1 = {1, 2, 3};
    const freeaml::Vector<int> v2 = {8, 26}; /* M * v1 */

    EXPECT_EQ(M * v1, v2);
}

TEST(SparseMatrixTest, RandomIntegerMatrixWithinDefaultRange)
{
    const freeaml::SparseMatrix<int> M =
        freeaml::random_sparse_matrix<int>(3, 4, 5);

    EXPECT_EQ(M.num_rows(), 3u);
    EXPECT_EQ(M.num_cols(), 4u);
    EXPECT_EQ(M.num_stored(), 5u);

    const freeaml::Vector<int> elements = M.flatten();

    for (const int x : elements)
    {
        EXPECT_GE(x, 0);
        EXPECT_LE(x, 1);
    }

    EXPECT_EQ(std::count(elements.begin(), elements.end(), 0), 7u);
}

TEST(SparseMatrixTest, RandomIntegerMatrixWithinSpecifiedRange)
{
    const freeaml::SparseMatrix<int> M =
        freeaml::random_sparse_matrix<int>(3, 4, 5, -5, 5);

    EXPECT_EQ(M.num_rows(), 3u);
    EXPECT_EQ(M.num_cols(), 4u);
    EXPECT_EQ(M.num_stored(), 5u);

    const freeaml::Vector<int> elements = M.flatten();

    for (const int x : elements)
    {
        EXPECT_GE(x, -5);
        EXPECT_LE(x, 5);
    }

    EXPECT_EQ(std::count(elements.begin(), elements.end(), 0), 7u);
}

TEST(SparseMatrixTest, RandomFloatingPointMatrixWithinDefaultRange)
{
    const freeaml::SparseMatrix<double> M =
        freeaml::random_sparse_matrix<double>(3, 4, 5);

    EXPECT_EQ(M.num_rows(), 3u);
    EXPECT_EQ(M.num_cols(), 4u);
    EXPECT_EQ(M.num_stored(), 5u);

    const freeaml::Vector<double> elements = M.flatten();

    for (const double x : elements)
    {
        EXPECT_GE(x, 0.0);
        EXPECT_LE(x, 1.0);
    }

    EXPECT_EQ(std::count(elements.begin(), elements.end(), 0.0), 7u);
}

TEST(SparseMatrixTest, RandomFloatingPointMatrixWithinSpecifiedRange)
{
    const freeaml::SparseMatrix<double> M =
        freeaml::random_sparse_matrix<double>(3, 4, 5, -5.0, 5.0);

    EXPECT_EQ(M.num_rows(), 3u);
    EXPECT_EQ(M.num_cols(), 4u);
    EXPECT_EQ(M.num_stored(), 5u);

    const freeaml::Vector<double> elements = M.flatten();

    for (const double x : elements)
    {
        EXPECT_GE(x, -5.0);
        EXPECT_LE(x, 5.0);
    }

    EXPECT_EQ(std::count(elements.begin(), elements.end(), 0.0), 7u);
}
