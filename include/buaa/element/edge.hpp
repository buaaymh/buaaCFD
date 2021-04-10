// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_EDGE_HPP_
#define INCLUDE_BUAA_ELEMENT_EDGE_HPP_

#include <cmath>

#include "buaa/element/node.hpp"

namespace buaa {
namespace element {

template <int kDim, class EdgeData, class CellData>
class Triangle;

template <class EdgeData, class CellData>
class Edge {
 public:
  // Types:
  using PointType = Point<2>;
  using NodeType = Node<2>;
  using CellType = Triangle<2, EdgeData, CellData>;
  using Data = EdgeData;
  // Constructors:
  Edge() = default;
  Edge(const NodeType& head, const NodeType& tail) :
      head_(head), tail_(tail) {}
  // Methods:
  const NodeType& Head() const { return head_; }
  const NodeType& Tail() const { return tail_; }
  Scalar Measure() const { return (head_ - tail_).norm(); }
  PointType Center() const { return PointType((head_ + tail_) * 0.5); }
  // Accessors:
  CellType* GetPositiveSide() const { return positive_side_; }
  CellType* GetNegativeSide() const { return negative_side_; }
  // Mutators:
  void SetPositiveSide(CellType* cell) { positive_side_ = cell; }
  void SetNegativeSide(CellType* cell) { negative_side_ = cell; }

  // Data::
  Data data;

 private:
  const NodeType& head_;
  const NodeType& tail_;
  CellType* positive_side_{nullptr};
  CellType* negative_side_{nullptr};
};


}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_EDGE_HPP_
