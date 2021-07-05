// Copyright 2021 Weicheng Pei and Minghao Yang

#include <iostream>

#include <Eigen/Dense>
#include "gtest/gtest.h"

namespace buaa {
namespace algebra {
using namespace Eigen;
class TestEigen : public ::testing::Test {
};
TEST_F(TestEigen, GetData) {
  Matrix3d a;
  a << 1, 2, 3, 1, 0, -1, 0, 1, 1;
  double* data = a.data();
  EXPECT_EQ(data[0], 1);
  EXPECT_EQ(data[1], 1);
  EXPECT_EQ(data[2], 0);
  EXPECT_EQ(data[3], 2);
  EXPECT_EQ(data[4], 0);
  EXPECT_EQ(data[5], 1);
  EXPECT_EQ(data[6], 3);
  EXPECT_EQ(data[7], -1);
  EXPECT_EQ(data[8], 1);
}
TEST_F(TestEigen, VectorArithmetic) {
  Vector3d a(0.0, 1.0, 2.0);
  EXPECT_EQ(a.rows(), 3);
  EXPECT_EQ(a.cols(), 1);
  Vector3d b;
  for (int i = 0; i < 3; ++i) {
    b(i) = i;
  }
  b += a;
  for (int i = 0; i < 3; ++i) {
    EXPECT_DOUBLE_EQ(b(i), i * 2);
  }
  b /= 2;
  for (int i = 0; i < 3; ++i) {
    EXPECT_DOUBLE_EQ(b(i), a(i));
  }
  EXPECT_DOUBLE_EQ(a.dot(b), 5);
}

TEST_F(TestEigen, MatrixArithmetic) {
  Matrix3d a;
  a << 1, 2, 3, 1, 0, -1, 0, 1, 1;
  auto b = a.inverse() * a;
  for (int i = 0; i < 3; ++i) {
    EXPECT_DOUBLE_EQ(b(i, i), 1);
  }
  Vector3d c(0, 1, 2);
  c = b * c;
  for (int i = 0; i < 3; ++i) {
    EXPECT_DOUBLE_EQ(c(i), i);
  }
}
}  // namespace algebra
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
