// Copyright 2021 Weicheng Pei and Minghao Yang
#include "buaa/element/triangle.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace element {

class TriangleTest : public ::testing::Test {
 protected:
  using T0 = element::Triangle<0>;
  using T1 = element::Triangle<1>;
  using T2 = element::Triangle<2>;
  using T3 = element::Triangle<3>;
  using NodeType = element::Node<2>;
  using PointType = element::Point<2>;
  NodeType a{0, 0.0, 0.0}, b{1, 1.0, 0.0}, c{2, 0.0, 2.0};
  Id id{1};
  Scalar eps{1e-6};
};
TEST_F(TriangleTest, ZeroDegreeCell) {
  auto cell = T0(id, a, b, c);
  EXPECT_EQ(cell.I(), id);
  EXPECT_EQ(cell.CountVertices(), 3);
  EXPECT_EQ(cell.A().X(), a.X());
  EXPECT_EQ(cell.GetNode(0).X(), a.X());
  EXPECT_EQ(cell.B().X(), b.X());
  EXPECT_EQ(cell.GetNode(1).X(), b.X());
  EXPECT_EQ(cell.C().X(), c.X());
  EXPECT_EQ(cell.GetNode(2).X(), c.X());
  EXPECT_EQ(cell.Measure(), 1.0);
  EXPECT_EQ(cell.Center().X() * 3, 1.0);
  EXPECT_EQ(cell.Center().Y() * 3, 2.0);
  EXPECT_EQ(cell.Degree(), 0);
}
TEST_F(TriangleTest, OneDegreeCell) {
  auto cell = T1(id, a, b, c);
  EXPECT_EQ(cell.I(), id);
  EXPECT_EQ(cell.Measure(), 1.0);
  EXPECT_EQ(cell.Degree(), 1);
  EXPECT_EQ(cell.DxInv(), 2.0);
  EXPECT_EQ(cell.DyInv(), 1.0);
  auto x_c = cell.Center().X();
  auto y_c = cell.Center().Y();
  // One Degree:
  Scalar x = 3.0, y = 2.0;
  EXPECT_EQ(cell.F_0_0_0(x, y), (x - x_c) * cell.DxInv());
  EXPECT_EQ(cell.F_0_1_0(x, y), 2);
  EXPECT_EQ(cell.F_1_0_0(x, y), (y - y_c) * cell.DyInv());
  EXPECT_EQ(cell.F_1_0_1(x, y), 1);
}
TEST_F(TriangleTest, TwoDegreeCell) {
  auto cell = T2(id, a, b, c);
  EXPECT_EQ(cell.I(), id);
  EXPECT_EQ(cell.Measure(), 1.0);
  EXPECT_EQ(cell.Degree(), 2);
  EXPECT_EQ(cell.DxInv(), 2.0);
  EXPECT_EQ(cell.DyInv(), 1.0);
  auto x_c = cell.Center().X();
  auto y_c = cell.Center().Y();
  // Two Degree:
  auto xx = 0;
  cell.Integrate([&](auto point){ return cell.F_2_0_0(point.X(), point.Y());}, &xx);
  EXPECT_NEAR(xx, 0, eps);
  EXPECT_EQ(cell.F_2_1_0(x_c, y_c), 0);
  EXPECT_EQ(cell.F_2_2_0(x_c, y_c), 8);
  auto xy = 0;
  cell.Integrate([&](auto point){ return cell.F_3_0_0(point.X(), point.Y());}, &xy);
  EXPECT_NEAR(xy, 0, eps);
  EXPECT_EQ(cell.F_3_1_0(x_c, y_c), 0);
  EXPECT_EQ(cell.F_3_0_1(x_c, y_c), 0);
  EXPECT_EQ(cell.F_3_1_1(x_c, y_c), 2);
  auto yy = 0;
  cell.Integrate([&](auto point){ return cell.F_4_0_0(point.X(), point.Y());}, &yy);
  EXPECT_NEAR(yy, 0, eps);
  EXPECT_EQ(cell.F_4_0_1(x_c, y_c), 0);
  EXPECT_EQ(cell.F_4_0_2(x_c, y_c), 2);
}
TEST_F(TriangleTest, ThreeDegreeCell) {
  auto cell = T3(id, a, b, c);
  EXPECT_EQ(cell.I(), id);
  EXPECT_EQ(cell.Measure(), 1.0);
  EXPECT_EQ(cell.Degree(), 3);
  EXPECT_EQ(cell.DxInv(), 2.0);
  EXPECT_EQ(cell.DyInv(), 1.0);
  auto x_c = cell.Center().X();
  auto y_c = cell.Center().Y();
  // Three Degree:
  auto xxx = 0;
  cell.Integrate([&](auto point){ return cell.F_5_0_0(point.X(), point.Y());}, &xxx);
  EXPECT_NEAR(xxx, 0, eps);
  EXPECT_EQ(cell.F_5_1_0(x_c, y_c), 0);
  EXPECT_EQ(cell.F_5_2_0(x_c, y_c), 0);
  EXPECT_EQ(cell.F_5_3_0(x_c, y_c), 48);
  auto xxy = 0;
  cell.Integrate([&](auto point){ return cell.F_6_0_0(point.X(), point.Y());}, &xxy);
  EXPECT_NEAR(xxy, 0, eps);
  EXPECT_EQ(cell.F_6_1_0(x_c, y_c), 0);
  EXPECT_EQ(cell.F_6_2_0(x_c, y_c), 0);
  EXPECT_EQ(cell.F_6_0_1(x_c, y_c), 0);
  EXPECT_EQ(cell.F_6_1_1(x_c, y_c), 0);
  EXPECT_EQ(cell.F_6_2_1(x_c, y_c), 8);
  auto xyy = 0;
  cell.Integrate([&](auto point){ return cell.F_5_0_0(point.X(), point.Y());}, &xyy);
  EXPECT_NEAR(xyy, 0, eps);
  EXPECT_EQ(cell.F_7_1_0(x_c, y_c), 0);
  EXPECT_EQ(cell.F_7_0_1(x_c, y_c), 0);
  EXPECT_EQ(cell.F_7_0_2(x_c, y_c), 0);
  EXPECT_EQ(cell.F_7_1_1(x_c, y_c), 0);
  EXPECT_EQ(cell.F_7_1_2(x_c, y_c), 4);
  auto yyy = 0;
  cell.Integrate([&](auto point){ return cell.F_8_0_0(point.X(), point.Y());}, &yyy);
  EXPECT_NEAR(yyy, 0, eps);
  EXPECT_EQ(cell.F_8_0_1(x_c, y_c), 0);
  EXPECT_EQ(cell.F_8_0_2(x_c, y_c), 0);
  EXPECT_EQ(cell.F_8_0_3(x_c, y_c), 6);
}
TEST_F(TriangleTest, LocalToGlobalXY) {
  auto cell = T3(id, a, b, c);
  auto global_p = cell.GetGlobalXY(1.0/3, 1.0/3, 1.0/3);
  EXPECT_EQ(global_p.X(), cell.Center().X());
  EXPECT_EQ(global_p.Y(), cell.Center().Y());
  auto global_a = cell.GetGlobalXY(1.0, 0.0, 0.0);
  EXPECT_EQ(global_a.X(), a.X());
  EXPECT_EQ(global_a.Y(), a.Y());
  auto global_b = cell.GetGlobalXY(0.0, 1.0, 0.0);
  EXPECT_EQ(global_b.X(), b.X());
  EXPECT_EQ(global_b.Y(), b.Y());
  auto global_c = cell.GetGlobalXY(0.0, 0.0, 1.0);
  EXPECT_EQ(global_c.X(), c.X());
  EXPECT_EQ(global_c.Y(), c.Y());
}
TEST_F(TriangleTest, Factorial) {
  auto cell = T3(id, a, b, c);
  EXPECT_EQ(cell.Factorial(0), 1);
  EXPECT_EQ(cell.Factorial(1), 1);
  EXPECT_EQ(cell.Factorial(2), 2);
  EXPECT_EQ(cell.Factorial(3), 6);
  EXPECT_EQ(cell.Factorial(4), 24);
}

}  // namespace element
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
