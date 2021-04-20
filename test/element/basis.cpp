// Copyright 2021 Minghao Yang
#include <iostream>

#include "buaa/element/basis.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace element {

class BasisFunctionTest : public ::testing::Test {
 protected:
  using T1 = element::Triangle<1>;
  using T2 = element::Triangle<2>;
  using T3 = element::Triangle<3>;
  using Node = element::Node<2>;
  using Point = element::Point<2>;
  Node a{0, 0.0, 0.0}, b{1, 1.0, 0.0}, c{2, 0.0, 2.0};
};
TEST_F(BasisFunctionTest, OneDegree) {
  T1 cell{1, a, b, c};
  using Basis = element::Basis<1>;
  auto basis = Basis(&cell);
  EXPECT_EQ(basis.Degree(), 1);
  EXPECT_EQ(basis.Center().X() * 3, 1);
  EXPECT_EQ(basis.Center().Y() * 3, 2);
  EXPECT_EQ(basis.Measure(), 1.0);
  EXPECT_EQ(basis.DxInv(), 2.0);
  EXPECT_EQ(basis.DyInv(), 1.0);
  auto x_c = basis.Center().X();
  auto y_c = basis.Center().Y();
  EXPECT_EQ(basis.F_0_0_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_0_1_0(x_c, y_c), 2);
  EXPECT_EQ(basis.F_1_0_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_1_0_1(x_c, y_c), 1);
}
TEST_F(BasisFunctionTest, TwoDegree) {
  T2 cell{1, a, b, c};
  using Basis = element::Basis<2>;
  auto basis = Basis(&cell);
  EXPECT_EQ(basis.Degree(), 2);
  EXPECT_EQ(basis.Center().X() * 3, 1);
  EXPECT_EQ(basis.Center().Y() * 3, 2);
  EXPECT_EQ(basis.Measure(), 1.0);
  EXPECT_EQ(basis.DxInv(), 2.0);
  EXPECT_EQ(basis.DyInv(), 1.0);
  auto x_c = basis.Center().X();
  auto y_c = basis.Center().Y();
  // One Degree:
  EXPECT_EQ(basis.F_0_0_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_0_1_0(x_c, y_c), 2);
  EXPECT_EQ(basis.F_1_0_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_1_0_1(x_c, y_c), 1);
  // Two Degree:
  auto xx = cell.Integrate([&](auto point){
    return std::pow(basis.F_0_0_0(point.X(), point.Y()), 2);}) / cell.Measure();
  EXPECT_EQ(basis.F_2_0_0(x_c, y_c), -xx);
  EXPECT_EQ(basis.F_2_1_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_2_2_0(x_c, y_c), 8);
  auto xy = cell.Integrate([&](auto point){
    return basis.F_0_0_0(point.X(), point.Y()) * basis.F_1_0_0(point.X(), point.Y());}) / cell.Measure();
  EXPECT_EQ(basis.F_3_0_0(x_c, y_c), -xy);
  EXPECT_EQ(basis.F_3_1_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_3_0_1(x_c, y_c), 0);
  EXPECT_EQ(basis.F_3_1_1(x_c, y_c), 2);
  auto yy = cell.Integrate([&](auto point){
    return std::pow(basis.F_1_0_0(point.X(), point.Y()), 2);}) / cell.Measure();
  EXPECT_EQ(basis.F_4_0_0(x_c, y_c), -yy);
  EXPECT_EQ(basis.F_4_0_1(x_c, y_c), 0);
  EXPECT_EQ(basis.F_4_0_2(x_c, y_c), 2);
}
TEST_F(BasisFunctionTest, ThreeDegree) {
  T3 cell{1, a, b, c};
  using Basis = element::Basis<3>;
  auto basis = Basis(&cell);
  EXPECT_EQ(basis.Degree(), 3);
  EXPECT_EQ(basis.Center().X() * 3, 1);
  EXPECT_EQ(basis.Center().Y() * 3, 2);
  EXPECT_EQ(basis.Measure(), 1.0);
  EXPECT_EQ(basis.DxInv(), 2.0);
  EXPECT_EQ(basis.DyInv(), 1.0);
  auto x_c = basis.Center().X();
  auto y_c = basis.Center().Y();
  // One Degree:
  EXPECT_EQ(basis.F_0_0_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_0_1_0(x_c, y_c), 2);
  EXPECT_EQ(basis.F_1_0_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_1_0_1(x_c, y_c), 1);
  // Two Degree:
  auto xx = cell.Integrate([&](auto point){
    return std::pow(basis.F_0_0_0(point.X(), point.Y()), 2);}) / cell.Measure();
  EXPECT_EQ(basis.F_2_0_0(x_c, y_c), -xx);
  EXPECT_EQ(basis.F_2_1_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_2_2_0(x_c, y_c), 8);
  auto xy = cell.Integrate([&](auto point){
    return basis.F_0_0_0(point.X(), point.Y()) * basis.F_1_0_0(point.X(), point.Y());}) / cell.Measure();
  EXPECT_EQ(basis.F_3_0_0(x_c, y_c), -xy);
  EXPECT_EQ(basis.F_3_1_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_3_0_1(x_c, y_c), 0);
  EXPECT_EQ(basis.F_3_1_1(x_c, y_c), 2);
  auto yy = cell.Integrate([&](auto point){
    return std::pow(basis.F_1_0_0(point.X(), point.Y()), 2);}) / cell.Measure();
  EXPECT_EQ(basis.F_4_0_0(x_c, y_c), -yy);
  EXPECT_EQ(basis.F_4_0_1(x_c, y_c), 0);
  EXPECT_EQ(basis.F_4_0_2(x_c, y_c), 2);
  // Three Degree:
  auto xxx = cell.Integrate([&](auto point){
    return std::pow(basis.F_0_0_0(point.X(), point.Y()), 3);}) / cell.Measure();
  EXPECT_EQ(basis.F_5_0_0(x_c, y_c), -xxx);
  EXPECT_EQ(basis.F_5_1_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_5_2_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_5_3_0(x_c, y_c), 48);
  auto xxy = cell.Integrate([&](auto point){
    return std::pow(basis.F_0_0_0(point.X(), point.Y()), 2) * basis.F_1_0_0(point.X(), point.Y()); }) / cell.Measure();
  EXPECT_EQ(basis.F_6_0_0(x_c, y_c), -xxy);
  EXPECT_EQ(basis.F_6_1_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_6_2_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_6_0_1(x_c, y_c), 0);
  EXPECT_EQ(basis.F_6_1_1(x_c, y_c), 0);
  EXPECT_EQ(basis.F_6_2_1(x_c, y_c), 8);
  auto xyy = cell.Integrate([&](auto point){
    return basis.F_0_0_0(point.X(), point.Y()) * std::pow(basis.F_1_0_0(point.X(), point.Y()), 2); }) / cell.Measure();
  EXPECT_EQ(basis.F_7_0_0(x_c, y_c), -xxy);
  EXPECT_EQ(basis.F_7_1_0(x_c, y_c), 0);
  EXPECT_EQ(basis.F_7_0_1(x_c, y_c), 0);
  EXPECT_EQ(basis.F_7_0_2(x_c, y_c), 0);
  EXPECT_EQ(basis.F_7_1_1(x_c, y_c), 0);
  EXPECT_EQ(basis.F_7_1_2(x_c, y_c), 4);
  auto yyy = cell.Integrate([&](auto point){
    return std::pow(basis.F_1_0_0(point.X(), point.Y()), 3);}) / cell.Measure();
  EXPECT_EQ(basis.F_8_0_0(x_c, y_c), -yyy);
  EXPECT_EQ(basis.F_8_0_1(x_c, y_c), 0);
  EXPECT_EQ(basis.F_8_0_2(x_c, y_c), 0);
  EXPECT_EQ(basis.F_8_0_3(x_c, y_c), 6);
}

}  // namespace element
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}