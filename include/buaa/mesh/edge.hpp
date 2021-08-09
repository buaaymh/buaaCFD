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
  static constexpr int nCoef = (kDegree+1) * (kDegree+2) / 2 - 1;
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // Types:
  using Cell = Triangle<kDegree, EdgeData, CellData>;
  using Base = element::Edge<kDegree>;
  using Point = element::Point<2>;
  using Node = element::Node<2>;
  using Matrix = typename Base::Matrix;
  using Data = EdgeData;
  // Constructors:
  Edge() = default;
  Edge(const Node& head, const Node& tail) : Base(head, tail) {
    b_matrix = Matrix::Zero();
  }
  // Accessors:
  Cell* GetPositiveSide() const { return positive_side_; }
  Cell* GetNegativeSide() const { return negative_side_; }
  Cell* GetOpposite(Cell* cell) const {
    if (cell == positive_side_) { return negative_side_; }
    else { return positive_side_; }
  }
  static constexpr int CountCoef() { return nCoef; }
  Scalar GetNormalX() const { return (Base::Tail().Y() - Base::Head().Y()) / Base::Measure(); }
  Scalar GetNormalY() const { return (Base::Head().X() - Base::Tail().X()) / Base::Measure(); }
  // Mutators:
  void SetPositiveSide(Cell* cell) { positive_side_ = cell; }
  void SetNegativeSide(Cell* cell) { negative_side_ = cell; }
  // Initialize VR Matrix and Vector:
  void InitializeBmat() {
    Scalar normal[2] = {GetNormalX(), GetNormalY()};
    this->Integrate([&](const Point& point) {
      return positive_side_->GetMatAt(point.X(), point.Y(), *negative_side_,
                                      distance, normal);
    }, &b_matrix);
  }
  void InitializeBmat(const Point& vec_ab) {
    Scalar normal[2] = {GetNormalX(), GetNormalY()};
    if (positive_side_->Contains(this)) {
      this->Integrate([&](const Point& point) {
        return positive_side_->GetMatAt(point.X(), point.Y(), *negative_side_,
                                                 distance, vec_ab, normal);
      }, &b_matrix);
    } else {
      this->Integrate([&](const Point& point) {
        return negative_side_->GetMatAt(point.X(), point.Y(), *positive_side_,
                                                 distance, vec_ab, normal);
      }, &b_matrix);
    }
  }
  // Data:
  Scalar distance;
  Data data;
  Matrix b_matrix;

 private:
  Cell* positive_side_{nullptr};
  Cell* negative_side_{nullptr};
};



}  // namespace mesh
}  // namespace buaa

#endif  //  INCLUDE_BUAA_MESH_EDGE_HPP_
