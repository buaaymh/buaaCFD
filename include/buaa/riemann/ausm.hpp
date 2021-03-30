// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_RIEMANN_AUSM_HPP_
#define INCLUDE_BUAA_RIEMANN_AUSM_HPP_

#include "buaa/riemann/types.hpp"

namespace buaa {
namespace riemann {

template <class GasModel, int kDim = 1>
class Ausm;

template <class GasModel>
class Ausm<GasModel, 2> {
 public:
  // Types:
  using Gas = GasModel;
  using FluxType = Flux<2>;
  using ConservativeType = Conservative<2>;
  using PrimitiveType = Primitive<2>;
  using State = PrimitiveType;
  using Vector = typename State::Vector;
  // Get F on T Axia
  FluxType GetFluxOnTimeAxis(State const& left, State const& right) {
    FluxType flux_positive = GetPositiveFlux(left);
    FluxType flux_negative = GetNegativeFlux(right);
    flux_positive += flux_negative;
    return flux_positive;
  }
  // Get F of U
  FluxType GetFlux(State const& state) {
    auto rho_u = state.rho() * state.u();
    auto rho_v = state.rho() * state.v();
    auto rho_u_u = rho_u * state.u();
    return {rho_u, rho_u_u + state.p(), rho_v * state.u(),
            state.u() * (state.p() * Gas::GammaOverGammaMinusOne()
                       + 0.5 * (rho_u_u + rho_v * state.v()))};
  }

 private:
  FluxType GetPositiveFlux(State const& state) {
    double p_positive = state.p();
    double a = Gas::GetSpeedOfSound(state);
    double mach = state.u() / a;
    double mach_positive = mach;
    double h = a * a / Gas::GammaMinusOne() + state.u() * state.u() * 0.5 +
                                              state.v() * state.v() * 0.5;
    FluxType flux = {1, state.u(), state.v(), h};
    if (mach >= -1 && mach <= 1) {
      mach_positive = (mach + 1) * (mach + 1) * 0.25;
      p_positive = state.p() * (mach + 1) * 0.5;
    } else if (mach < -1) {
      mach_positive = 0.0;
      p_positive = 0.0;
    }
    double temp = state.rho() * a * mach_positive;
    flux *= temp;
    flux.momentum(0) += p_positive;
    return flux;
  }
  FluxType GetNegativeFlux(State state) {
    double p_negative = state.p();
    double a = Gas::GetSpeedOfSound(state);
    double mach = state.u() / a;
    double mach_negative = mach;
    double h = a * a / Gas::GammaMinusOne() + state.u() * state.u() * 0.5 +
                                              state.v() * state.v() * 0.5;
    FluxType flux = {1, state.u(), state.v(), h};
    if (mach >= -1 && mach <= 1) {
      mach_negative = - (mach - 1) * (mach - 1) * 0.25;
      p_negative = - state.p() * (mach - 1) * 0.5;
    } else if (mach > 1) {
      mach_negative = 0.0;
      p_negative = 0.0;
    }
    double temp = state.rho() * a * mach_negative;
    flux *= temp;
    flux.momentum(0) += p_negative;
    return flux;
  }
};

}  //  namespace riemann
}  //  namespace buaa

#endif  //  INCLUDE_BUAA_RIEMANN_AUSM_HPP_
