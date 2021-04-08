// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_MESH_EDGE_HPP_
#define INCLUDE_BUAA_MESH_EDGE_HPP_

#include <cmath>

#include "buaa/element/edge.hpp"
#include "buaa/mesh/triangle.hpp"

namespace buaa {
namespace mesh {

using Scalar = double;

template <int kDim, class EdgeData, class CellData>
class Triangle;

template <class EdgeData, class CellData>
class Edge : public element::Edge {
 public:
  // Types:
  using Node = element::Edge::NodeType;
  using Cell = mesh::Triangle<2, EdgeData, CellData>;
  using Data = EdgeData;
  // Constructors:
  template <class... Args>
  explicit Edge(Args&&... args) :
      element::Edge(std::forward<Args>(args)...) {}
  Data data;

 private:
  Cell* positive_side_{nullptr};
  Cell* negative_side_{nullptr};
};

}  // namespace mesh
}  // namespace buaa

#endif  //  INCLUDE_BUAA_MESH_EDGE_HPP_
