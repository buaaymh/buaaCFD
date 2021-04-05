// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_NODE_HPP_
#define INCLUDE_BUAA_ELEMENT_NODE_HPP_

#include <array>
#include <cstddef>

namespace buaa {
namespace element {
using Scalar = double;
using Id = std::size_t;

template <int kDim>
class Node;

template <>
class Node<2> {
 public:
  // Constructors:
  Node() = default;
  Node(Scalar x, Scalar y) : coordinates_({x, y}) {}
  // Methods:
  Scalar X() const { return coordinates_[0]; }
  Scalar Y() const { return coordinates_[1]; }

 private:
  std::array<Scalar, 2> coordinates_;
};

template <>
class Node<3> {
 public:
  // Constructors:
  Node() = default;
  Node(Scalar x, Scalar y, Scalar z) : coordinates_({x, y, z}) {}

  Scalar X() const { return coordinates_[0]; }
  Scalar Y() const { return coordinates_[1]; }
  Scalar Z() const { return coordinates_[2]; }

 private:
  std::array<Scalar, 3> coordinates_;
};

}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_NODE_HPP_
