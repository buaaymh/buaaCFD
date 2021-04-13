// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_
#define INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_

#include <cmath>
#include <initializer_list>
#include <vector>

#include <Eigen/Dense>

#include "buaa/element/edge.hpp"
#include "buaa/element/node.hpp"

namespace buaa {
namespace element {

template <class EdgeData, class CellData>
class Edge;

template <int kDim, class EdgeData, class CellData>
class Triangle {
 public:
  // Types:
  using Edge = element::Edge<EdgeData, CellData>;
  using Node = typename Edge::NodeType;
  using Point = typename Edge::PointType;;
  using Data = CellData;
  // Constructors:
  Triangle() = default;
  Triangle(const Node& a, const Node& b, const Node& c,
      std::initializer_list<Edge*> edges) : a_(a), b_(b), c_(c), edges_{edges} {
      auto cross = (b_.X() - a_.X()) * (c_.Y() - a_.Y()) -
                   (b_.Y() - a_.Y()) * (c_.X() - a_.X());   
      measure_ = std::abs(cross) * 0.5;
      center_ = Point((a + b + c) / 3);
  }
  // Accessors:
  int CountVertices() const { return 3; }
  const Node& A() const { return a_; }
  const Node& B() const { return b_; }
  const Node& C() const { return c_; }
  const Node& GetNode(int i) const {
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
  const Point& Center() const { return center_; }
  // Data:
  Data data;

 private:
  const Node& a_;
  const Node& b_;
  const Node& c_;
  std::vector<Edge*> edges_;
  Scalar measure_;
  Point center_;
};

}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_