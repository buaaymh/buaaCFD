// Copyright 2021 Weicheng Pei and Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_DATA_HPP_
#define INCLUDE_BUAA_ELEMENT_DATA_HPP_

#include <array>

#include "buaa/element/point.hpp"

namespace buaa {
namespace mesh {

using Scalar = element::Scalar;
using Id = std::size_t;

template <int kDim, int kScalars, int kVectors>
struct Data {
 public:
  using Vector = std::array<Scalar, kDim>;
  static constexpr int CountScalars() { return kScalars; }
  static constexpr int CountVectors() { return kVectors; }
  std::array<Scalar, kScalars> scalars;
  std::array<Vector, kVectors> vectors;
};

using Empty = Data<0, 0, 0>;

}  // namespace mesh
}  // namespace buaa

#endif  // INCLUDE_BUAA_ELEMENT_DATA_HPP_
