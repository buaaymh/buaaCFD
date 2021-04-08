// Copyright 2021 Minghao Yang

#include "buaa/mesh/data.hpp"
#include "buaa/mesh/edge.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace mesh {

class EdgeTest : public ::testing::Test {
 protected:
  struct EdgeData {
   public:
    Scalar measure;
    std::array<Scalar, 2> center;
  };
  using Edge = mesh::Edge<EdgeData, Empty>;
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
TEST_F(EdgeTest, Data) {
  auto edge = Edge(head, tail);
  edge.data.measure = edge.Measure();
  edge.data.center = {edge.Center().X(), edge.Center().Y()};
  EXPECT_EQ(edge.data.measure, edge.Measure());
  auto center = edge.Center();
  EXPECT_EQ(edge.data.center[0] * 2, head.X() + tail.X());
  EXPECT_EQ(edge.data.center[1] * 2, head.Y() + tail.Y());
}

}  // namespace mesh
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
