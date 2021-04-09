// Copyright 2021 Weicheng Pei and Minghao Yang
#include "buaa/element/point.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace element {

class PointTest : public ::testing::Test {
 protected:
  using P2 = Point<2>;
  using P3 = Point<3>;
  const Scalar x{1.0}, y{2.0}, z{3.0};
};
TEST_F(PointTest, TemplateConstructor) {
  // Test P2(Real x, Real y):
  auto p2 = P2(x, y);
  EXPECT_EQ(p2.X(), x);
  EXPECT_EQ(p2.Y(), y);
  // Test P3(Real x, Real y, Real z):
  auto p3 = P3(x, y, z);
  EXPECT_EQ(p3.X(), x);
  EXPECT_EQ(p3.Y(), y);
  EXPECT_EQ(p3.Z(), z);
}
TEST_F(PointTest, EqualOperator) {
  auto a = P3(x, y, z);
  auto b = P3(x, y, z);
  EXPECT_EQ(a, b);
  auto c = P3(0.0, 0.0, 0.0);
  EXPECT_NE(a, c);
}
TEST_F(PointTest, PlusAndMinusOperator) {
  auto a = P3(x, y, z);
  auto b = P3(0.0, 0.0, 0.0);
  EXPECT_EQ(a + b, a);
  EXPECT_EQ(a - b, a);
  EXPECT_EQ(b + b, b);
  EXPECT_EQ(a - a, b);
}
TEST_F(PointTest, ScalarMultiplication) {
  auto a = P3(x, y, z);
  auto b = P3(0.0, 0.0, 0.0);
  EXPECT_EQ(a * 1, a);
  EXPECT_EQ(a * 0, b);
  EXPECT_EQ(a / 1, a);
  EXPECT_EQ(b / 2, b);
}
TEST_F(PointTest, DotProduct) {
  auto a = P3(x, y, z);
  auto b = P3(0.0, 0.0, 0.0);
  EXPECT_EQ(a.dot(a), 14);
  EXPECT_EQ(a.dot(b), 0);
}

}  // namespace element
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
