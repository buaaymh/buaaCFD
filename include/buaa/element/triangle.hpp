// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_
#define INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_

#include <cmath>
#include <stdexcept>

#include <Eigen/Dense>

#include "buaa/element/node.hpp"

namespace buaa {
namespace element {

template <int kDim>
class Triangle;

template <>
class Triangle<2> {
 public:
  // Types:
  using NodeType = Node<2>;
  // Constructors:
  Triangle() = default;
  Triangle(const NodeType& a, const NodeType& b, const NodeType& c)
      : a_(a), b_(b), c_(c) {
      auto cross = (b_.X() - a_.X()) * (c_.Y() - a_.Y()) -
                   (b_.Y() - a_.Y()) * (c_.X() - a_.X());
      measure_ = std::abs(cross) * 0.5;
      center_ = NodeType((a_.X() + b_.X() + c_.X()) / 3,
                         (a_.Y() + b_.Y() + c_.Y()) / 3);
  }
  // Accessors:
  int CountVertices() const { return 3; }
  const NodeType& A() const { return a_; }
  const NodeType& B() const { return b_; }
  const NodeType& C() const { return c_; }
  const NodeType& GetPoint(int i) const {
    switch (i)  {
    case 0:
      return A();
    case 1:
      return B();
    case 2:
      return C();
    default:
      throw std::out_of_range("A `Triangle` has 3 `Point`s.");
    }
  }
  // Geometric methods:
  Scalar Measure() const {return measure_; }
  const NodeType& Center() const { return center_; }

 private:
  const NodeType& a_;
  const NodeType& b_;
  const NodeType& c_;
  Scalar measure_;
  NodeType center_;
};

}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_