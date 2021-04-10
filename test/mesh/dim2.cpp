// Copyright 2021 Minghao Yang

#include <vector>

#include "buaa/mesh/dim2.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace mesh {

class MeshTest : public ::testing::Test {
 protected:
  using MeshType = Mesh<>;
  using CellType = MeshType::Cell;
  using EdgeType = MeshType::Edge;
  using NodeType = MeshType::Node;
  MeshType mesh{};
  const std::vector<Scalar> x{0.0, 1.0, 1.0, 0.0}, y{0.0, 0.0, 1.0, 1.0};
};
TEST_F(MeshTest, DefaultConstructor) {
  EXPECT_EQ(mesh.CountNodes(), 0);
  EXPECT_EQ(mesh.CountEdges(), 0);
  EXPECT_EQ(mesh.CountCells(), 0);
}
TEST_F(MeshTest, EmplaceNode) {
  mesh.SetNodesNum(1);
  EXPECT_EQ(mesh.CountNodes(), 0);
  mesh.EmplaceNode(0, 0.0, 0.0);
  EXPECT_EQ(mesh.CountNodes(), 1);
}
TEST_F(MeshTest, ForEachNode) {
  // Emplace 4 nodes:
  mesh.SetNodesNum(4);
  for (auto i = 0; i != x.size(); ++i) {
    mesh.EmplaceNode(i, x[i], y[i]);
  }
  EXPECT_EQ(mesh.CountNodes(), x.size());
  // Check each node's index and coordinates:
  mesh.ForEachNode([&](NodeType const& node) {
    auto i = node.I();
    EXPECT_EQ(node.X(), x[i]);
    EXPECT_EQ(node.Y(), y[i]);
  });
}
TEST_F(MeshTest, EmplaceEdge) {
  mesh.SetNodesNum(2);
  mesh.SetEdgesNum(1);
  EXPECT_EQ(mesh.CountNodes(), 0);
  mesh.EmplaceNode(0, 0.0, 0.0);
  mesh.EmplaceNode(1, 1.0, 0.0);
  EXPECT_EQ(mesh.CountNodes(), 2);
  EXPECT_EQ(mesh.CountEdges(), 0);
  mesh.EmplaceEdge(0, 1);
  EXPECT_EQ(mesh.CountEdges(), 1);
}
TEST_F(MeshTest, ForEachEdge) {
  /*
     3 ----- 2
     | \   / |
     |   X   |
     | /   \ |
     0 ----- 1
  */
  // Emplace 4 nodes:
  mesh.SetNodesNum(4);
  for (auto n = 0; n != x.size(); ++n) {
    mesh.EmplaceNode(n, x[n], y[n]);
  }
  EXPECT_EQ(mesh.CountNodes(), x.size());
  // Emplace 6 edges:
  mesh.SetEdgesNum(6);
  EXPECT_EQ(mesh.CountEdges(), 0);
  mesh.EmplaceEdge(0, 1);
  mesh.EmplaceEdge(1, 2);
  mesh.EmplaceEdge(2, 3);
  mesh.EmplaceEdge(3, 0);
  mesh.EmplaceEdge(2, 0);
  mesh.EmplaceEdge(3, 1);
  EXPECT_EQ(mesh.CountEdges(), 6);
  // For each edge: head's index < tail's index
  mesh.ForEachEdge([](EdgeType const& edge) {
    EXPECT_LT(edge.Head().I(), edge.Tail().I());
  });
}
TEST_F(MeshTest, EmplaceCell) {
  /*
     3 ----- 2
     | (0) / |
     |   /   |
     | / (1) |
     0 ----- 1
  */
  // Emplace 4 nodes:
  for (auto n = 0; n != x.size(); ++n) {
    mesh.EmplaceNode(n, x[n], y[n]);
  }
  EXPECT_EQ(mesh.CountNodes(), x.size());
  // Emplace 2 triangular cells:
  mesh.SetEdgesNum(5);
  mesh.SetCellsNum(2);
  EXPECT_EQ(mesh.CountEdges(), 0);
  EXPECT_EQ(mesh.CountCells(), 0);
  mesh.EmplaceCell(0, {0, 1, 2});
  mesh.EmplaceCell(1, {0, 2, 3});
  EXPECT_EQ(mesh.CountCells(), 2);
  EXPECT_EQ(mesh.CountEdges(), 5);
}
TEST_F(MeshTest, GetSide) {
  /*
     3 -- [2] -- 2
     |  (1)   /  |
    [3]   [4]   [1]
     |  /   (0)  |
     0 -- [0] -- 1
  */
  mesh.Clear();
  // Emplace 4 nodes:
  for (auto n = 0; n != x.size(); ++n) {
    mesh.EmplaceNode(n, x[n], y[n]);
  }
  EXPECT_EQ(mesh.CountNodes(), x.size());
  // Emplace 5 edges:
  auto edges = std::vector<EdgeType*>();
  edges.emplace_back(mesh.EmplaceEdge(0, 1));
  edges.emplace_back(mesh.EmplaceEdge(1, 2));
  edges.emplace_back(mesh.EmplaceEdge(2, 3));
  edges.emplace_back(mesh.EmplaceEdge(3, 0));
  edges.emplace_back(mesh.EmplaceEdge(0, 2));
  EXPECT_EQ(mesh.CountEdges(), edges.size());
  // Emplace 2 triangular cells:
  auto cells = std::vector<CellType*>();
  cells.emplace_back(mesh.EmplaceCell(0, {0, 1, 2}));
  cells.emplace_back(mesh.EmplaceCell(1, {0, 2, 3}));
  EXPECT_EQ(mesh.CountCells(), 2);
  EXPECT_EQ(mesh.CountEdges(), edges.size());
  // Check each edge's positive side and negative side:
  // edges[0] == {nodes[0], nodes[1]}
  EXPECT_EQ(edges[0]->GetPositiveSide(), cells[0]);
  EXPECT_EQ(edges[0]->GetNegativeSide(), nullptr);
  // edges[1] == {nodes[1], nodes[2]}
  EXPECT_EQ(edges[1]->GetPositiveSide(), cells[0]);
  EXPECT_EQ(edges[1]->GetNegativeSide(), nullptr);
  // edges[4] == {nodes[0], nodes[2]}
  EXPECT_EQ(edges[4]->GetPositiveSide(), cells[1]);
  EXPECT_EQ(edges[4]->GetNegativeSide(), cells[0]);
  // edges[2] == {nodes[2], nodes[3]}
  EXPECT_EQ(edges[2]->GetPositiveSide(), cells[1]);
  EXPECT_EQ(edges[2]->GetNegativeSide(), nullptr);
  // edges[3] == {nodes[0], nodes[3]}
  EXPECT_EQ(edges[3]->GetPositiveSide(), nullptr);
  EXPECT_EQ(edges[3]->GetNegativeSide(), cells[1]);
}


}  // namespace mesh
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}