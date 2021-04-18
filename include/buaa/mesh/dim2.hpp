// Copyright 2021 Minghao Yang

#ifndef INCLUDE_BUAA_MESH_DIM2_HPP_
#define INCLUDE_BUAA_MESH_DIM2_HPP_

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <forward_list>
#include <initializer_list>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "buaa/element/node.hpp"
#include "buaa/element/edge.hpp"
#include "buaa/element/triangle.hpp"
#include "buaa/element/data.hpp"

namespace buaa {
namespace mesh {

using Id = element::Id;
using Scalar = element::Scalar;
using Empty = element::Empty;

template <class EdgeData = Empty,
          class CellData = Empty>
class Mesh;

template <class EdgeData, class CellData>
class Mesh {
 public:
  using Cell = typename element::Triangle<2, EdgeData, CellData>;
  using Edge = typename Cell::Edge;
  using Node = typename Cell::Node;
  using Point = typename Cell::Point;
  using NodeId = Id;
  using EdgeId = Id;
  using CellId = Id;

  // Constructors:
  Mesh() = default;
  void SetNodesNum(int num) { id_to_node_.reserve(num); }
  void SetEdgesNum(int num) { id_to_edge_.reserve(num); }
  void SetCellsNum(int num) { id_to_cell_.reserve(num); }
  // Count primitive objects.
  auto CountNodes() const { return id_to_node_.size(); }
  auto CountEdges() const { return id_to_edge_.size(); }
  auto CountCells() const { return id_to_cell_.size(); }
  // Traverse primitive objects.
  template <class Visitor>
  void ForEachNode(Visitor&& visitor) const {
    for (auto& node_ptr : id_to_node_) { visitor(*node_ptr); }
  }
  template <class Visitor>
  void ForEachEdge(Visitor&& visitor) const {
    for (auto& edge_ptr : id_to_edge_) { visitor(*edge_ptr); }
  }
  template <class Visitor>
  void ForEachCell(Visitor&& visitor) const {
    for (auto& cell_ptr : id_to_cell_) { visitor(*cell_ptr); }
  }
  // Emplace primitive objects.
  Node* EmplaceNode(NodeId i, Scalar x, Scalar y) {
    auto node_unique_ptr = std::make_unique<Node>(i, x, y);
    auto node_ptr = node_unique_ptr.get();
    id_to_node_.emplace(id_to_node_.begin()+i, std::move(node_unique_ptr));
    return node_ptr;
  }
  Edge* EmplaceEdge(NodeId head_id, NodeId tail_id) {
    if (head_id > tail_id) { std::swap(head_id, tail_id); }
    auto node_pair = std::minmax(head_id, tail_id);
    auto iter = node_pair_to_edge_.find(node_pair);
    if (iter != node_pair_to_edge_.end()) {
      return iter->second;
    } else {  // Emplace a new edge:
      auto edge_unique_ptr = std::make_unique<Edge>(*(id_to_node_.at(head_id)),
                                                    *(id_to_node_.at(tail_id)));
      auto edge_ptr = edge_unique_ptr.get();
      node_pair_to_edge_.emplace(node_pair, edge_ptr);
      id_to_edge_.emplace_back(std::move(edge_unique_ptr));
      assert(id_to_edge_.size() == node_pair_to_edge_.size());
      return edge_ptr;
    }
  }
  Node* GetNode(NodeId i) const { return id_to_node_.at(i).get(); }
  Cell* EmplaceCell(CellId i, std::initializer_list<NodeId> nodes) {
    auto* p = nodes.begin();
    auto* a = GetNode(p[0]);
    auto* b = GetNode(p[1]);
    auto* c = GetNode(p[2]);
    if (IsClockWise(*a, *b, *c)) {
      std::swap(a, c);
    }
    auto ab = EmplaceEdge(p[0], p[1]);
    auto bc = EmplaceEdge(p[1], p[2]);
    auto ca = EmplaceEdge(p[2], p[0]);
    auto edges = {ab, bc, ca};
    auto cell_unique_ptr = std::make_unique<Cell>(i, *a, *b, *c, edges);
    auto cell_ptr = cell_unique_ptr.get();
    id_to_cell_.emplace_back(std::move(cell_unique_ptr));
    LinkCellToEdge(cell_ptr, p[0], p[1], ab);
    LinkCellToEdge(cell_ptr, p[1], p[2], bc);
    LinkCellToEdge(cell_ptr, p[2], p[0], ca);
    return cell_ptr;
  }
  void Clear() {
    id_to_node_.clear();
    id_to_edge_.clear();
    id_to_cell_.clear();
    node_pair_to_edge_.clear();
  }
  static constexpr int Dim() { return 2; }

 private:
  void LinkCellToEdge(Cell* cell, NodeId head_id, NodeId tail_id,
                      Edge* edge) {
    if (head_id < tail_id) {
      edge->SetPositiveSide(cell);
    } else {
      edge->SetNegativeSide(cell);
    }
  }
  bool IsClockWise(Node const& a, Node const& b, Node const& c) {
    auto ab = Point(b - a);
    auto ac = Point(c - a);
    auto cross = ab.X() * ac.Y() - ab.Y() * ac.X();
    return cross < 0;
  }

 private:
  std::vector<std::unique_ptr<Node>> id_to_node_;
  std::vector<std::unique_ptr<Edge>> id_to_edge_;
  std::vector<std::unique_ptr<Cell>> id_to_cell_;
  std::map<std::pair<NodeId, NodeId>, Edge*> node_pair_to_edge_;
};

}  // namespace mesh
}  // namespace buaa

#endif  // INCLUDE_BUAA_MESH_DIM2_HPP_