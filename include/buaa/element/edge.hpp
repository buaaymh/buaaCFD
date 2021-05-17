// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_EDGE_HPP_
#define INCLUDE_BUAA_ELEMENT_EDGE_HPP_

#include <cmath>

#include "buaa/riemann/types.hpp"
#include "buaa/element/gauss.hpp"
#include "buaa/element/node.hpp"

namespace buaa {
namespace element {

template <int kDegree>
class Edge {
 private:
  static constexpr int num_coefficients = (kDegree+1) * (kDegree+2) / 2 - 1;
  static constexpr int num_quad_points = int(kDegree / 2) + 1;
 public:
  // Types:
  using PointType = Point<2>;
  using NodeType = Node<2>;
  using Gauss = element::Gauss<num_quad_points>;
  using Matrix = Eigen::Matrix<Scalar, num_coefficients, num_coefficients>;
  using Vector = Eigen::Matrix<Scalar, num_coefficients, 1>;
  using Flux = riemann::Flux<2>;
  // Constructors:
  Edge() = default;
  Edge(const NodeType& head, const NodeType& tail) : head_(head), tail_(tail) {
    for (int i = 0; i < num_quad_points; ++i) {
      quad_points_[i] = LocalToGlobal(gauss_.x_local[i]);
    }
  }
  // Accessors:
  static constexpr int Degree() { return kDegree; }
  static constexpr int CountQuadPoints() { return num_quad_points; }
  static Gauss GetGauss() { return gauss_; }
  // Methods:
  const NodeType& Head() const { return head_; }
  const NodeType& Tail() const { return tail_; }
  Scalar Measure() const { return (head_ - tail_).norm(); }
  PointType Center() const { return PointType((head_ + tail_) * 0.5); }
  template <class Visitor>
  void ForEachQuadPoint(Visitor&& visitor) const {
    for (auto& point : quad_points_) { visitor(point); }
  }
  template <class Integrand>
  Scalar Integrate(Integrand&& integrand) const {
    Scalar result{0.0};
    for (int i = 0; i < num_quad_points; ++i) {
      result += integrand(quad_points_[i]) * gauss_.weights[i];
    }
    result *= 0.5 * Measure();
    return result;
  }
  // template <class Value, class Integrand>
  // Value Integrate(Integrand&& integrand) const {
  //   Scalar result;
  //   for (int i = 0; i < num_quad_points; ++i) {
  //     result += integrand(quad_points_[i]) * gauss_.weights[i];
  //   }
  //   result *= 0.5 * Measure();
  //   return result;
  // }
  template <class Integrand>
  Matrix IntegrateM(Integrand&& integrand) const {
    Matrix result;
    for (int i = 0; i < num_quad_points; ++i) {
      result += integrand(quad_points_[i]) * gauss_.weights[i];
    }
    return result * 0.5 * Measure();
  }
  template <class Integrand>
  Vector IntegrateV(Integrand&& integrand) const {
    Vector result;
    for (int i = 0; i < num_quad_points; ++i) {
      result += integrand(quad_points_[i]) * gauss_.weights[i];
    }
    return result * 0.5 * Measure();
  }

 private:
  PointType LocalToGlobal(Scalar x_local) {
    auto point = ((Head() + Tail()) + (Tail() - Head()) * x_local) * 0.5;
    return PointType{point(0), point(1)};
  }
 private:
  const NodeType& head_;
  const NodeType& tail_;
  std::array<PointType, num_quad_points> quad_points_;
  static constexpr Gauss gauss_ = Gauss();
};

}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_EDGE_HPP_
