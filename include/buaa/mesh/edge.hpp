// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_MESH_EDGE_HPP_
#define INCLUDE_BUAA_MESH_EDGE_HPP_

#include <array>
#include <utility>

#include <Eigen/Dense>

#include "buaa/mesh/data.hpp"
#include "buaa/element/edge.hpp"
// #include "buaa/mesh/dim2.hpp"
#include "buaa/mesh/gauss.hpp"
// #include "buaa/mesh/triangle.hpp"

namespace buaa {
namespace mesh {

template <int kDegree, class CellData, class EdgeData>
class Triangle;

template <int   kDegree,
          class EdgeData = Empty,
          class CellData = Empty>
class Edge;

template <int kDegree, class EdgeData, class CellData>
class Edge : public element::Edge {
 private:
  static constexpr int num_quad_points = int(kDegree / 2) + 1;
  static constexpr int num_coefficients = (kDegree+1) * (kDegree+2) / 2 - 1;
 public:
  // Types:
  using Base = element::Edge;
  using Point = element::Point<2>;
  using Node = element::Node<2>;
  using Cell = Triangle<kDegree, CellData, EdgeData>;
  using Gauss = mesh::Gauss<num_quad_points>;
  using Matrix = Eigen::Matrix<Scalar, num_coefficients, num_coefficients>;
  using Data = EdgeData;
  // Constructors:
  Edge() = default;
  Edge(const Node& head, const Node& tail) : Base(head, tail) {
    for (int i = 0; i < num_quad_points; ++i) {
      quad_points_[i] = LocalToGlobal(gauss_.x_local[i]);
    }
  }
  // Accessors:
  Cell* GetPositiveSide() const { return positive_side_; }
  Cell* GetNegativeSide() const { return negative_side_; }
  // // Mutators:
  void SetPositiveSide(Cell* cell) { positive_side_ = cell; }
  void SetNegativeSide(Cell* cell) { negative_side_ = cell; }
  // Accessors:
  static constexpr int Degree() { return kDegree; }
  static constexpr int CountQuadPoints() { return num_quad_points; }
  static constexpr int CountCoefficients() { return num_coefficients; }
  static Gauss GetGauss() { return gauss_; }
  // Methods:
  template <class Visitor>
  void ForEachQuadPoint(Visitor&& visitor) const {
    for (auto& point : quad_points_) { visitor(point); }
  }
  template <class Integrand>
  Scalar Integrate(Integrand&& integrand) const {
    Scalar result = 0.0;
    for (int i = 0; i < num_quad_points; ++i) {
      result += integrand(quad_points_[i]) * gauss_.weights[i];
    }
    return result * 0.5 * Measure();
  }
  // Data:
  Data data;
  Matrix b_matrix;

 private:
  Point LocalToGlobal(Scalar x_local) {
    auto point = ((Head() + Tail()) + (Tail() - Head()) * x_local) * 0.5;
    return Point{point(0), point(1)};
  }

 private:
  Cell* positive_side_{nullptr};
  Cell* negative_side_{nullptr};
  static constexpr Gauss gauss_ = Gauss();
  std::array<Point, num_quad_points> quad_points_;
};

}  // namespace mesh
}  // namespace buaa

#endif  //  INCLUDE_BUAA_MESH_EDGE_HPP_
