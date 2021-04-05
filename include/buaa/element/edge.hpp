// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_EDGE_HPP_
#define INCLUDE_BUAA_ELEMENT_EDGE_HPP_

#include <cmath>

#include "buaa/element/node.hpp"

namespace buaa {
namespace element {

class Edge {
 public:
  // Types:
  using NodeType = Node<2>;
  // Constructors:
  Edge() = default;
  Edge(const NodeType& head, const NodeType& tail) :
      head_(head), tail_(tail) {}
  // Methods:
  const NodeType& Head() const { return head_; }
  const NodeType& Tail() const { return tail_; }
  Scalar Measure() const {
    auto dx = head_.X() - tail_.X();
    auto dy = head_.Y() - tail_.Y();
    return std::sqrt(dx * dx + dy * dy);
  }
  NodeType Center() const {
    auto x_c = (head_.X() + tail_.X()) * 0.5;
    auto y_c = (head_.Y() + tail_.Y()) * 0.5;
    return NodeType{x_c, y_c};
  }

 private:
  const NodeType& head_;
  const NodeType& tail_;
};


}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_EDGE_HPP_
