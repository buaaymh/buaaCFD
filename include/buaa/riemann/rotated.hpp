// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_RIEMANN_ROTATED_HPP_
#define INCLUDE_BUAA_RIEMANN_ROTATED_HPP_

#include "buaa/riemann/ausm.hpp"
#include "buaa/riemann/linear.hpp"
#include "buaa/riemann/types.hpp"

namespace buaa {
namespace riemann {

template <class UnrotatedSimple>
class Simple {
  using Base = UnrotatedSimple;

 public:
  using Vector = typename Base::Vector;
  using State = typename Base::State;
  using Flux = typename Base::Flux;
  using Coefficient = Vector;
  void Rotate(Vector const& normal) {
    Rotate(normal(0), normal(1));
  }
  void Rotate(Scalar const& n_1, Scalar const& n_2) {
    auto a_normal = global_coefficient(0) * n_1;
    a_normal += global_coefficient(1) * n_2;
    unrotated_simple_ = UnrotatedSimple(a_normal);
  }
  Flux GetFluxOnTimeAxis(State const& left, State const& right) {
    auto flux = unrotated_simple_.GetFluxOnTimeAxis(left, right);
    return flux;
  }
  Flux GetFluxOnSolidWall(State const& state) {
    return {};
  }
  Flux GetFluxOnFreeWall(State const& state) {
    return unrotated_simple_.GetFlux(state);
  }
  static Coefficient global_coefficient;

 private:
  UnrotatedSimple unrotated_simple_;
};
template <class UnrotatedSimple>
typename Simple<UnrotatedSimple>::Coefficient
Simple<UnrotatedSimple>::global_coefficient;

template <class UnrotatedEuler>
class RotatedEuler {
  using Base = UnrotatedEuler;

 public:
  using Gas = typename Base::Gas;
  using Vector = typename Base::Vector;
  using FluxType = typename Base::FluxType;
  using ConservativeType = typename Base::ConservativeType;
  using PrimitiveType = typename Base::PrimitiveType;
  using State = ConservativeType;
  void Rotate(Vector const& normal) { normal_ = normal; }
  void Rotate(Scalar const& n_1, Scalar const& n_2) {
    normal_(0) = n_1;
    normal_(1) = n_2;
  }
  FluxType GetFluxOnTimeAxis(
      ConservativeType const& left,
      ConservativeType const& right) {
    auto left__primitive = Gas::ConservativeToPrimitive(left);
    auto right_primitive = Gas::ConservativeToPrimitive(right);
    GlobalToNormal(&(left__primitive.momentum));
    GlobalToNormal(&(right_primitive.momentum));
    auto flux = unrotated_euler_.GetFluxOnTimeAxis(
        left__primitive, right_primitive);
    NormalToGlobal(&(flux.momentum));
    return flux;
  }
  FluxType GetFluxOnSolidWall(ConservativeType const& conservative) {
    auto primitive = Gas::ConservativeToPrimitive(conservative);
    auto flux = FluxType();
    flux.momentum(0) = primitive.p();
    NormalToGlobal(&(flux.momentum));
    return flux;
  }
  FluxType GetFluxOnFreeWall(ConservativeType const& conservative) {
    auto primitive = Gas::ConservativeToPrimitive(conservative);
    GlobalToNormal(&(primitive.momentum));
    auto flux = unrotated_euler_.GetFlux(primitive);
    NormalToGlobal(&(flux.momentum));
    return flux;
  }
  void GlobalToNormal(Vector* v) {
    auto& n = normal_;
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
  void NormalToGlobal(Vector* v) {
    auto& n = normal_;
    auto v_0 = (*v)(0) * n(0) - (*v)(1) * n(1);
    (*v)(1) = (*v)(0) * n(1) + (*v)(1) * n(0);
    (*v)(0) = v_0;
  }

 private:
  UnrotatedEuler unrotated_euler_;
  Vector normal_;
};

}  // namespace riemann
}  // namespace buaa

#endif  //  INCLUDE_BUAA_RIEMANN_ROTATED_HPP_
