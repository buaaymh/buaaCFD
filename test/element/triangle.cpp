// Copyright 2021 Weicheng Pei and Minghao Yang
#include "buaa/element/triangle.hpp"
#include "buaa/element/data.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace element {

class TriangleTest : public ::testing::Test {
 protected:
  using TriangleType = Triangle;
  using EdgeType = TriangleType::EdgeType;
  using NodeType = TriangleType::NodeType;
  using PointType = TriangleType::PointType;
  NodeType a{0, 0.0, 0.0}, b{1, 3.0, 0.0}, c{2, 0.0, 3.0};
  EdgeType ab{a, b}, bc{b, c}, ca{c, a};
  Id id{1};
};
TEST_F(TriangleTest, TemplateConstructor) {
  auto cell = TriangleType(id, a, b, c, {&ab, &bc, &ca});
  EXPECT_EQ(cell.I(), id);
  EXPECT_EQ(cell.CountVertices(), 3);
}
TEST_F(TriangleTest, GeometricMethods) {
  auto cell = TriangleType(id, a, b, c, {&ab, &bc, &ca});
  EXPECT_EQ(cell.I(), id);
  EXPECT_EQ(cell.A().X(), a.X());
  EXPECT_EQ(cell.GetNode(0).X(), a.X());
  EXPECT_EQ(cell.B().X(), b.X());
  EXPECT_EQ(cell.GetNode(1).X(), b.X());
  EXPECT_EQ(cell.C().X(), c.X());
  EXPECT_EQ(cell.GetNode(2).X(), c.X());
  EXPECT_EQ(cell.Measure(), 4.5);
  EXPECT_EQ(cell.Center().X(), 1.0);
  EXPECT_EQ(cell.Center().Y(), 1.0);
}
TEST_F(TriangleTest, LocalToGlobalXY) {
  auto cell = TriangleType(id, a, b, c, {&ab, &bc, &ca});
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
TEST_F(TriangleTest, Integrator) {
  auto cell = TriangleType(id, a, b, c, {&ab, &bc, &ca});
  auto integrand_0 = [](auto x, auto y) { return 3.14; };
  EXPECT_DOUBLE_EQ(cell.Integrate(integrand_0), cell.Measure() * 3.14);
  auto integrand_1 = [](auto x, auto y) { return x - 1.0; };
  EXPECT_DOUBLE_EQ(cell.Integrate(integrand_1), 0.0);
  auto integrand_2 = [](auto x, auto y) { return y - 1.0; };
  EXPECT_DOUBLE_EQ(cell.Integrate(integrand_2), 0.0);
  auto integrand_3 = [](auto x, auto y) { return x * y; };
  EXPECT_DOUBLE_EQ(cell.Integrate(integrand_3), 3.375);
  auto integrand_4 = [](auto x, auto y) { return x * x * y; };
  EXPECT_DOUBLE_EQ(cell.Integrate(integrand_4), 4.05);
}

}  // namespace element
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
