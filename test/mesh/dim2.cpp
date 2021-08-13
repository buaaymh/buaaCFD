// Copyright 2021 Minghao Yang

#include <vector>

#include "buaa/mesh/dim2.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace mesh {

class MeshTest : public ::testing::Test {
 protected:
  using MeshType = Mesh<1, Empty, Empty>;
  using CellType = MeshType::Cell;
  using EdgeType = MeshType::Edge;
  using PointType = MeshType::Point;
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
TEST_F(MeshTest, ForEachCell) {
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
  // For each cell: head's index < tail's index
  Id id{0};
  mesh.ForEachCell([&](CellType const& cell) {
    EXPECT_EQ(cell.I(), id++);
    EXPECT_EQ(cell.Measure(), 0.5);
  });
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
  // Add Ghost Cells
  auto move = PointType{1, 0};
  auto ghost_l = *cells[0]; ghost_l.Move(move);
  auto ghost_r = *cells[1]; ghost_l.Move(move);
  // Check each edge's positive side and negative side:
  // edges[0] == {nodes[0], nodes[1]}
  EXPECT_EQ(edges[0]->GetPositiveSide(), cells[0]);
  EXPECT_EQ(edges[0]->GetNegativeSide(), nullptr);
  EXPECT_EQ(edges[0]->GetOpposite(nullptr), cells[0]);
  EXPECT_EQ(edges[0]->GetOpposite(cells[0]), nullptr);
  // edges[1] == {nodes[1], nodes[2]}
  EXPECT_EQ(edges[1]->GetPositiveSide(), cells[0]);
  EXPECT_EQ(edges[1]->GetNegativeSide(), nullptr);
  EXPECT_EQ(edges[1]->GetOpposite(nullptr), cells[0]);
  EXPECT_EQ(edges[1]->GetOpposite(cells[0]), nullptr);
  // edges[4] == {nodes[0], nodes[2]}
  EXPECT_EQ(edges[4]->GetPositiveSide(), cells[1]);
  EXPECT_EQ(edges[4]->GetNegativeSide(), cells[0]);
  EXPECT_EQ(edges[4]->GetOpposite(cells[0]), cells[1]);
  EXPECT_EQ(edges[4]->GetOpposite(cells[1]), cells[0]);
  // edges[2] == {nodes[2], nodes[3]}
  EXPECT_EQ(edges[2]->GetPositiveSide(), cells[1]);
  EXPECT_EQ(edges[2]->GetNegativeSide(), nullptr);
  EXPECT_EQ(edges[2]->GetOpposite(nullptr), cells[1]);
  EXPECT_EQ(edges[2]->GetOpposite(cells[1]), nullptr);
  // edges[3] == {nodes[0], nodes[3]}
  EXPECT_EQ(edges[3]->GetPositiveSide(), nullptr);
  EXPECT_EQ(edges[3]->GetNegativeSide(), cells[1]);
  EXPECT_EQ(edges[3]->GetOpposite(nullptr), cells[1]);
  EXPECT_EQ(edges[3]->GetOpposite(cells[1]), nullptr);
}
TEST_F(MeshTest, GetMatrix) {
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
  // Emplace 5 edges:
  auto edges = std::vector<EdgeType*>();
  edges.emplace_back(mesh.EmplaceEdge(0, 1));
  edges.emplace_back(mesh.EmplaceEdge(1, 2));
  edges.emplace_back(mesh.EmplaceEdge(2, 3));
  edges.emplace_back(mesh.EmplaceEdge(3, 0));
  edges.emplace_back(mesh.EmplaceEdge(0, 2));
  // Emplace 2 triangular cells:
  auto cells = std::vector<CellType*>();
  cells.emplace_back(mesh.EmplaceCell(0, {0, 1, 2}));
  cells.emplace_back(mesh.EmplaceCell(1, {0, 2, 3}));
  EXPECT_EQ(mesh.CountNodes(), x.size());
  EXPECT_EQ(mesh.CountCells(), 2);
  EXPECT_EQ(mesh.CountEdges(), 5);
  auto x_31 = PointType(edges[1]->Center() - edges[3]->Center());
  auto x_02 = PointType(edges[2]->Center() - edges[0]->Center());
  edges[0]->distance = (edges[0]->GetPositiveSide()->Center() -
                        cells[1]->Center() + x_02).norm();
  edges[1]->distance = (edges[1]->GetPositiveSide()->Center() -
                        cells[1]->Center() - x_31).norm();
  edges[2]->distance = (edges[2]->GetPositiveSide()->Center() -
                        cells[0]->Center() - x_02).norm();
  edges[3]->distance = (cells[0]->Center() - x_31 -
                        edges[3]->GetNegativeSide()->Center()).norm();
  edges[4]->distance = (edges[4]->GetPositiveSide()->Center() -
                        edges[4]->GetNegativeSide()->Center()).norm();
  edges[4]->InitializeBmat();
  auto b_mat = edges[4]->b_matrix;
  EXPECT_EQ(b_mat.size(), 4);
  std::cout << b_mat << std::endl << std::endl;
  cells[0]->InitializeAmatInv();
  auto a_mat_inv = cells[0]->a_matrix_inv;
  EXPECT_EQ(a_mat_inv.size(), 4);
  std::cout << a_mat_inv << std::endl << std::endl;
  cells[1]->InitializeAmatInv();
  a_mat_inv = cells[1]->a_matrix_inv;
  EXPECT_EQ(a_mat_inv.size(), 4);
  std::cout << a_mat_inv << std::endl << std::endl;
}

}  // namespace mesh
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}