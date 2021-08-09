// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_GAUSS_HPP_
#define INCLUDE_BUAA_ELEMENT_GAUSS_HPP_

#include "buaa/element/point.hpp"
#include <array>

namespace buaa {
namespace element {

template <int kPoint>
struct Gauss;

template <>
struct Gauss<1> {
 public:
  Gauss() = default;
  int CountPoint() const { return 1; }
  std::array<Scalar, 1> x_local{0.0};
  std::array<Scalar, 1> weights{2.0};
};
template <>
struct Gauss<2> {
 public:
  Gauss() = default;
  int CountPoint() const { return 2; }
  std::array<Scalar, 2> x_local{-0.5773502691896250, 0.5773502691896250};
  std::array<Scalar, 2> weights{1.0, 1.0};
};

}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_GAUSS_HPP_