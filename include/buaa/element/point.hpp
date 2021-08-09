// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_POINT_HPP_
#define INCLUDE_BUAA_ELEMENT_POINT_HPP_

#include <Eigen/Dense>

#include "buaa/riemann/types.hpp"

namespace buaa {
namespace element {
using Scalar = riemann::Scalar;
using Id = std::size_t;
using Eigen::Matrix;

template <int kDim>
class Point : public Matrix<Scalar, kDim, 1> {
 public:
  // Types:
  using Base = Matrix<Scalar, kDim, 1>;
  // Constructors:
  Point() = default;
  template<class... XYZ>
  explicit Point(XYZ&&... xyz) : Base{std::forward<XYZ>(xyz)...} {}
  // Methods:
  Scalar X() const { return (*this)(0); }
  Scalar Y() const { return (*this)(1); }
  Scalar Z() const { return (*this)(2); }
};

}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_POINT_HPP_
