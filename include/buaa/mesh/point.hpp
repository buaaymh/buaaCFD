// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_MESH_POINT_HPP_
#define INCLUDE_BUAA_MESH_POINT_HPP_

#include <array>
#include <cstddef>

namespace buaa {
namespace mesh {
using Scalar = double;
using Id = std::size_t;

template <int kDim>
class Point;

template <>
class Point<2> {
 public:
  // Constructors:
  Point() = default;
  Point(Id id, Scalar x, Scalar y) : id_(id),
      coordinates_({x, y}) {}

  Id I() const { return id_; }
  Scalar X() const { return coordinates_[0]; }
  Scalar Y() const { return coordinates_[1]; }

 private:
  Id id_;
  std::array<Scalar, 2> coordinates_;
};

template <>
class Point<3> {
 public:
  // Constructors:
  Point() = default;
  Point(Id id, Scalar x, Scalar y, Scalar z) : id_(id),
      coordinates_({x, y, z}) {}

  Id I() const { return id_; }
  Scalar X() const { return coordinates_[0]; }
  Scalar Y() const { return coordinates_[1]; }
  Scalar Z() const { return coordinates_[2]; }

 private:
  Id id_;
  std::array<Scalar, 3> coordinates_;
};

}  // namespace mesh
}  // namespace buaa

#endif  //  INCLUDE_BUAA_MESH_POINT_HPP_
