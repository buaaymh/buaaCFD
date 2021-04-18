// Copyright 2021 Minghao Yang
#include "buaa/mesh/edge.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace mesh {

class EdgeTest : public ::testing::Test {
 protected:
  using E1 = mesh::Edge<1, Empty, Empty>;
  using E2 = mesh::Edge<2, Empty, Empty>;
  using E3 = mesh::Edge<3, Empty, Empty>;
  using Node = element::Node<2>;
  Node head{0, 0.3, 0.0}, tail{1, 0.0, 0.4};
};
TEST_F(EdgeTest, Constructor) {
  // Test Edge(Id, const Point &, const Point &):
  auto edge = E1(head, tail);
  EXPECT_EQ(edge.Head().X(), head.X());
  EXPECT_EQ(edge.Tail().X(), tail.X());
  EXPECT_EQ(edge.Head().Y(), head.Y());
  EXPECT_EQ(edge.Tail().Y(), tail.Y());
  EXPECT_EQ(edge.CountQuadPoints(), 1);
}
TEST_F(EdgeTest, QuadPoints) {
  auto e1 = E1(head, tail);
  EXPECT_EQ(e1.GetOrder(), 1);
  EXPECT_EQ(e1.CountQuadPoints(), 1);
  EXPECT_EQ(e1.CountCoefficients(), 2);
  auto e2 = E2(head, tail);
  EXPECT_EQ(e2.GetOrder(), 2);
  EXPECT_EQ(e2.CountQuadPoints(), 2);
  auto e3 = E3(head, tail);
  EXPECT_EQ(e3.GetOrder(), 3);
  EXPECT_EQ(e3.CountQuadPoints(), 2);
}

}  // namespace mesh
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}