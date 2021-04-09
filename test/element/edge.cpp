// Copyright 2021 Minghao Yang
#include "buaa/element/edge.hpp"
#include "buaa/mesh/data.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace element {

class EdgeTest : public ::testing::Test {
 protected:
  struct EdgeData {
   public:
    Scalar measure;
    std::array<Scalar, 2> center;
  };
  using Edge = element::Edge<EdgeData, mesh::Empty>;
  using Node = typename Edge::NodeType;
  Node head{0, 0.3, 0.0}, tail{1, 0.0, 0.4};
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
TEST_F(EdgeTest, Data) {
  auto edge = Edge(head, tail);
  edge.data.measure = edge.Measure();
  edge.data.center = {edge.Center().X(), edge.Center().Y()};
  EXPECT_EQ(edge.data.measure, edge.Measure());
  auto center = edge.Center();
  EXPECT_EQ(edge.data.center[0] * 2, head.X() + tail.X());
  EXPECT_EQ(edge.data.center[1] * 2, head.Y() + tail.Y());
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
