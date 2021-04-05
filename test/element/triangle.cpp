// Copyright 2021 Weicheng Pei and Minghao Yang
#include "buaa/element/triangle.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace element {

class TriangleTest : public ::testing::Test {
 protected:
  using NodeType = Node<2>;
  using TriangleType = Triangle<2>;
  NodeType a{0.0, 0.0}, b{1.0, 0.0}, c{0.0, 1.0};
};
TEST_F(TriangleTest, TemplateConstructor) {
  auto cell = TriangleType(a, b, c);
  EXPECT_EQ(cell.CountVertices(), 3);
}
TEST_F(TriangleTest, GeometricMethods) {
  auto cell = TriangleType(a, b, c);
  EXPECT_EQ(cell.A().X(), a.X());
  EXPECT_EQ(cell.GetPoint(0).X(), a.X());
  EXPECT_EQ(cell.B().X(), b.X());
  EXPECT_EQ(cell.GetPoint(1).X(), b.X());
  EXPECT_EQ(cell.C().X(), c.X());
  EXPECT_EQ(cell.GetPoint(2).X(), c.X());
  EXPECT_EQ(cell.Measure(), 0.5);
  auto node = cell.Center();
  EXPECT_EQ(node.X() * 3, a.X() + b.X() + c.X());
  EXPECT_EQ(node.Y() * 3, a.Y() + b.Y() + c.Y());
}

}  // namespace element
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
