// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_MESH_EDGE_HPP_
#define INCLUDE_BUAA_MESH_EDGE_HPP_

#include <array>
#include <utility>

#include <Eigen/Dense>

#include "buaa/element/edge.hpp"
#include "buaa/mesh/dim2.hpp"
// #include "buaa/mesh/triangle.hpp"

namespace buaa {
namespace mesh {

template <int Order, class CellData, class EdgeData>
class Triangle;

template <int   Order,
          class EdgeData = Empty,
          class CellData = Empty>
class Edge;

template <int Order, class EdgeData, class CellData>
class Edge : public element::Edge<EdgeData, CellData> {
 private:
  static constexpr int num_coefficients = (Order+1) * (Order+2) / 2 - 1;
  static constexpr int num_quad_points = int(Order/2)+1;
  
 public:
  // Types:
  using Base = element::Edge<EdgeData, CellData>;
  using Point = element::Point<2>;
  using Node = element::Node<2>;
  using Cell = Triangle<2, Order, CellData, EdgeData>;
  using Matrix = Eigen::Matrix<Scalar, num_coefficients, num_coefficients>;
  // Constructors:
  Edge() = default;
  Edge(const Node& head, const Node& tail) : Base(head, tail) {}
  // template<class... Args>
  // explicit Edge(Args&&... args) : Base{std::forward<Args>(args)...} {}

  // Accessors:
  Cell* GetPositiveSide() const { return positive_side_; }
  Cell* GetNegativeSide() const { return negative_side_; }
  // Mutators:
  void SetPositiveSide(Cell* cell) { positive_side_ = cell; }
  void SetNegativeSide(Cell* cell) { negative_side_ = cell; }
  // Accessors:
  static constexpr int GetOrder() { return Order; }
  static constexpr int CountQuadPoints() { return num_quad_points; }
  static constexpr int CountCoefficients() { return num_coefficients; }
  const Matrix& B() const { return b_matrix; }
  // Methods:
  template <class Visitor>
  void ForEachQuadPoint(Visitor&& visitor) const {
    for (auto& point : quad_points_) { visitor(point); }
  }

 private:
  std::array<Point, num_quad_points> quad_points_;
  Matrix b_matrix;
  Cell* positive_side_{nullptr};
  Cell* negative_side_{nullptr};

};


}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_MESH_EDGE_HPP_
