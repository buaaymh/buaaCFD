// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_RIEMANN_LINEAR_HPP_
#define INCLUDE_BUAA_RIEMANN_LINEAR_HPP_

#include "buaa/riemann/types.hpp"

namespace buaa {
namespace riemann {

class Linear {
 public:
  static constexpr int kDim = 2;
  // Types:
  using Vector = Matrix<Scalar, kDim, 1>;
  using Coefficient = Matrix<Scalar, kDim, 1>;
  using Flux = Scalar;
  using Jacobi = double;
  using State = Scalar;
  // Constructor:
  Linear() : a_const_(1) {}
  explicit Linear(Scalar const& a_const) : a_const_(a_const) {}
  // Get F on T Axia:
  Flux GetFluxOnTimeAxis(State const& left, State const& right) const {
    if (0 < a_const_) {
      return left * a_const_;
    } else {
      return right* a_const_;
    }
  }
  // Get F of U
  Flux GetFlux(State const& state) const { return state * a_const_ ; }

 private:
  Scalar a_const_;
};

}  // namespace riemann
}  // namespace buaa

#endif  // INCLUDE_BUAA_RIEMANN_LINEAR_HPP_
