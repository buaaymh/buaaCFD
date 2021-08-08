// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_RIEMANN_LINEAR_HPP_
#define INCLUDE_BUAA_RIEMANN_LINEAR_HPP_

#include "buaa/riemann/types.hpp"

namespace buaa {
namespace riemann {

class Linear {
 public:
  // Get F of U_l and U_r
  static Scalar GetFlux(Scalar const& left, Scalar const& right,
                           Scalar const& a) {
    if (0 < a) { return left * a; }
    else { return right* a; }
  }
  // Get F of U
  static Scalar GetFlux(Scalar const& state, Scalar const& a) {
    return state * a;
  }
};

}  // namespace riemann
}  // namespace buaa

#endif  // INCLUDE_BUAA_RIEMANN_LINEAR_HPP_
