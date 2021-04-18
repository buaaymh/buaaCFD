// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_EDGE_HPP_
#define INCLUDE_BUAA_ELEMENT_EDGE_HPP_

#include <cmath>

#include "buaa/element/node.hpp"
#include "buaa/element/data.hpp"

namespace buaa {
namespace element {

class Edge {
 public:
  // Types:
  using PointType = Point<2>;
  using NodeType = Node<2>;
  // Constructors:
  Edge() = default;
  Edge(const NodeType& head, const NodeType& tail) :
      head_(head), tail_(tail) {}
  // Methods:
  const NodeType& Head() const { return head_; }
  const NodeType& Tail() const { return tail_; }
  Scalar Measure() const { return (head_ - tail_).norm(); }
  PointType Center() const { return PointType((head_ + tail_) * 0.5); }

 private:
  const NodeType& head_;
  const NodeType& tail_;
};


}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_EDGE_HPP_
