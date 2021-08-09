// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_NODE_HPP_
#define INCLUDE_BUAA_ELEMENT_NODE_HPP_

#include <cmath>

#include "buaa/element/point.hpp"

namespace buaa {
namespace element {

using Id = std::size_t;

template <int kDim>
class Node : public Point<kDim> {
 public:
  // Constructors:
  template<class... XYZ>
  explicit Node(Id id, XYZ&&... xyz)
      : id_(id), Point<kDim>{std::forward<XYZ>(xyz)...} {}
  Id I() const { return id_; }
 private:
  Id id_;
};

}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_NODE_HPP_
