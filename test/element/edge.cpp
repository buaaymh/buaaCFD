// Copyright 2021 Minghao Yang
#include "buaa/element/edge.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace element {

class EdgeTest : public ::testing::Test {
 protected:
  using Node = typename Edge::NodeType;
  Node head{0.3, 0.0}, tail{0.0, 0.4};
};
TEST_F(EdgeTest, Constructor) {
  // Test Edge(Id, const Point &, const Point &):
  auto edge = Edge(head, tail);
  EXPECT_EQ(edge.Head().X(), head.X());
  EXPECT_EQ(edge.Tail().X(), tail.X());
  EXPECT_EQ(edge.Head().Y(), head.Y());
  EXPECT_EQ(edge.Tail().Y(), tail.Y());
}
TEST_F(EdgeTest, GeometryMethods) {
  auto edge = Edge(head, tail);
  EXPECT_EQ(edge.Measure(), 0.5);
  auto center = edge.Center();
  EXPECT_EQ(center.X() * 2, head.X() + tail.X());
  EXPECT_EQ(center.Y() * 2, head.Y() + tail.Y());
}
// TEST_F(EdgeTest, ElementMethods) {
//   auto edge = Edge(i, head, tail);
//   auto integrand = [](const auto& point) { return 3.14; };
//   EXPECT_EQ(edge.Integrate(integrand), edge.Measure() * 3.14);
// }

}  // namespace element
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
