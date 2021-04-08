// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_MESH_TRIANGLE_HPP_
#define INCLUDE_BUAA_MESH_TRIANGLE_HPP_

#include <array>
#include <cmath>
#include <initializer_list>
#include <vector>

#include "buaa/mesh/edge.hpp"
#include "buaa/element/triangle.hpp"

namespace buaa {
namespace mesh {

template <class EdgeData, class CellData>
class Edge;

template <int kDim, class EdgeData, class CellData>
class Triangle : public element::Triangle<kDim> {
 public:
  // Types:
  using Edge = mesh::Edge<EdgeData, CellData>;
  using Node = typename Edge::Node;
  using Data = CellData;
  // Constructors:
  Triangle(const Node& a, const Node& b, const Node& c,
      std::initializer_list<Edge*> edges) : element::Triangle<kDim>(a, b, c),
      edges_{edges} {}
  Data data;
 protected:
  std::vector<Edge*> edges_;
};

}  // namespace mesh
}  // namespace buaa

#endif  //  INCLUDE_BUAA_MESH_TRIANGLE_HPP_
