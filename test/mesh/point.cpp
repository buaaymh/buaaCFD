// Copyright 2021 Weicheng Pei and Minghao Yang
#include "buaa/mesh/point.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace mesh {

class PointTest : public ::testing::Test {
 protected:
  using P2 = Point<2>;
  using P3 = Point<3>;
  const int i{8};
  const Scalar x{1.0}, y{2.0}, z{3.0};
};
TEST_F(PointTest, TemplateConstructor) {
  // Test P2(Id i, Real x, Real y):
  auto p2 = P2(i, x, y);
  EXPECT_EQ(p2.I(), i);
  EXPECT_EQ(p2.X(), x);
  EXPECT_EQ(p2.Y(), y);
  // Test P3(Id i, Real x, Real y, Real z):
  auto p3 = P3(i, x, y, z);
  EXPECT_EQ(p3.I(), i);
  EXPECT_EQ(p3.X(), x);
  EXPECT_EQ(p3.Y(), y);
  EXPECT_EQ(p3.Z(), z);
}

}  // namespace mesh
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
