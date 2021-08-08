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
  static FluxType GetFlux(State const& left, State const& right) {
    FluxType flux_positive = GetPositiveFlux(left);
    FluxType flux_negative = GetNegativeFlux(right);
    flux_positive += flux_negative;
    return flux_positive;
  }
  // Get F on normal T Axia
  static FluxType GetFlux(ConservativeType const& left,
                          ConservativeType const& right,
                          Vector const& normal) {
    auto left__primitive = Gas::ConservativeToPrimitive(left);
    auto right_primitive = Gas::ConservativeToPrimitive(right);
    GlobalToNormal(&(left__primitive.momentum), normal);
    GlobalToNormal(&(right_primitive.momentum), normal);
    auto flux = GetFluxOnTimeAxis(left__primitive, right_primitive);
    NormalToGlobal(&(flux.momentum), normal);
    return flux;
  }
  // Get F of U
  static FluxType GetFlux(State const& state) {
    auto rho_u = state.rho() * state.u();
    auto rho_v = state.rho() * state.v();
    auto rho_u_u = rho_u * state.u();
    return {rho_u, rho_u_u + state.p(), rho_v * state.u(),
            state.u() * (state.p() * Gas::GammaOverGammaMinusOne()
                       + 0.5 * (rho_u_u + rho_v * state.v()))};
  }
  static void GlobalToNormal(Vector* v, Vector const& n) {
    /* Calculate the normal component: */
    auto v_n = v->dot(n);
    /* Calculate the tangential component:
       auto t = Vector{ -n[1], n[0] };
       auto v_t = v->Dot(t);
    */
    (*v)(1) = n(0) * (*v)(1) - n(1) * (*v)(0);
    /* Write the normal component: */
    (*v)(0) = v_n;
  }
  static void NormalToGlobal(Vector* v, Vector const& n) {
    auto v_0 = (*v)(0) * n(0) - (*v)(1) * n(1);
    (*v)(1) = (*v)(0) * n(1) + (*v)(1) * n(0);
    (*v)(0) = v_0;
  }
  
 private:
  static FluxType GetPositiveFlux(State const& state) {
    Scalar p_positive = state.p();
    Scalar a = Gas::GetSpeedOfSound(state);
    Scalar mach = state.u() / a;
    Scalar mach_positive = mach;
    Scalar h = a * a / Gas::GammaMinusOne() + (state.u() * state.u() +
                                               state.v() * state.v()) * 0.5;
    FluxType flux = {1, state.u(), state.v(), h};
    if (mach >= -1 && mach <= 1) {
      mach_positive = (mach + 1) * (mach + 1) * 0.25;
      p_positive = state.p() * (mach + 1) * 0.5;
    } else if (mach < -1) {
      mach_positive = 0.0;
      p_positive = 0.0;
    }
    Scalar temp = state.rho() * a * mach_positive;
    flux *= temp;
    flux.momentum(0) += p_positive;
    return flux;
  }
  static FluxType GetNegativeFlux(State const& state) {
    Scalar p_negative = state.p();
    Scalar a = Gas::GetSpeedOfSound(state);
    Scalar mach = state.u() / a;
    Scalar mach_negative = mach;
    Scalar h = a * a / Gas::GammaMinusOne() + (state.u() * state.u() +
                                               state.v() * state.v()) * 0.5;
    FluxType flux = {1, state.u(), state.v(), h};
    if (mach >= -1 && mach <= 1) {
      mach_negative = - (mach - 1) * (mach - 1) * 0.25;
      p_negative = - state.p() * (mach - 1) * 0.5;
    } else if (mach > 1) {
      mach_negative = 0.0;
      p_negative = 0.0;
    }
    Scalar temp = state.rho() * a * mach_negative;
    flux *= temp;
    flux.momentum(0) += p_negative;
    return flux;
  }
};

}  //  namespace riemann
}  //  namespace buaa

#endif  //  INCLUDE_BUAA_RIEMANN_AUSM_HPP_
