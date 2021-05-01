// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_MESH_EDGE_HPP_
#define INCLUDE_BUAA_MESH_EDGE_HPP_

#include <array>
#include <utility>

#include <Eigen/Dense>

#include "buaa/mesh/data.hpp"
#include "buaa/element/edge.hpp"
#include "buaa/mesh/triangle.hpp"

namespace buaa {
namespace mesh {

template <int kDegree, class EdgeData, class CellData>
class Triangle;

template <int   kDegree,
          class EdgeData = Empty,
          class CellData = Empty>
class Edge;

template <int kDegree, class EdgeData, class CellData>
class Edge : public element::Edge<kDegree> {
 private:
  static constexpr int num_coefficients = (kDegree+1) * (kDegree+2) / 2 - 1;
 public:
  // Types:
  using Cell = Triangle<kDegree, EdgeData, CellData>;
  using Base = element::Edge<kDegree>;
  using Point = element::Point<2>;
  using Node = element::Node<2>;
  using Matrix = Eigen::Matrix<Scalar, num_coefficients, num_coefficients>;
  using Data = EdgeData;
  // Constructors:
  Edge() = default;
  Edge(const Node& head, const Node& tail) : Base(head, tail) {}
  // Accessors:
  Cell* GetPositiveSide() const { return positive_side_; }
  Cell* GetNegativeSide() const { return negative_side_; }
  Cell* GetOpposite(Cell* cell) const {
    if (cell == positive_side_) { return negative_side_; }
    else { return positive_side_; }
  }
  // // Mutators:
  void SetPositiveSide(Cell* cell) { positive_side_ = cell; }
  void SetNegativeSide(Cell* cell) { negative_side_ = cell; }
  // Accessors:
  static constexpr int CountCoefficients() { return num_coefficients; }
  // Data:
  Data data;
  Matrix b_matrix;

 private:
  Cell* positive_side_{nullptr};
  Cell* negative_side_{nullptr};
};



}  // namespace mesh
}  // namespace buaa

#endif  //  INCLUDE_BUAA_MESH_EDGE_HPP_
