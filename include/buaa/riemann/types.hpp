// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_RIEMANN_TYPES_HPP_
#define INCLUDE_BUAA_RIEMANN_TYPES_HPP_

#include <cmath>

#include <Eigen/Dense>

namespace buaa {
namespace riemann {
using Scalar = double;
using Eigen::Matrix;

template <int kDim>
class Tuple {
 public:
  // Types:
  using Vector = Matrix<Scalar, kDim, 1>;
  // Data:
  Scalar mass{0};
  Vector momentum;
  Scalar energy{0};
  // Constructors:
  Tuple() = default;
  Tuple(Scalar const& rho,
        Scalar const& u,
        Scalar const& p)
      : mass{rho}, energy{p}, momentum(u) {}
  Tuple(Scalar const& rho,
        Scalar const& u, Scalar const& v,
        Scalar const& p)
      : mass{rho}, energy{p}, momentum(u, v) {}
  Tuple(Scalar const& rho,
        Scalar const& u, Scalar const& v, Scalar const& w,
        Scalar const& p)
      : mass{rho}, energy{p}, momentum(u, v, w) {}
  // Arithmetic Operators:
  Tuple& operator+=(Tuple const& that) {
    this->mass += that.mass;
    this->energy += that.energy;
    this->momentum += that.momentum;
    return *this;
  }
  Tuple& operator-=(Tuple const& that) {
    this->mass -= that.mass;
    this->energy -= that.energy;
    this->momentum -= that.momentum;
    return *this;
  }
  Tuple& operator*=(Scalar const& s) {
    this->mass *= s;
    this->energy *= s;
    this->momentum *= s;
    return *this;
  }
  Tuple& operator/=(Scalar const& s) {
    this->mass /= s;
    this->energy /= s;
    this->momentum /= s;
    return *this;
  }
  // Other Operations:
  bool operator==(Tuple const& that) const {
    return (this->mass == that.mass) && (this->energy == that.energy) &&
           (this->momentum == that.momentum);
  }
  Tuple operator*(const Scalar& s) {
    Tuple value = *this;
    value *= s;
    return value;
  }
};
template <int kDim>
class Flux : public Tuple<kDim> {
  // Types:
  using Base = Tuple<kDim>;
  // Constructors:
  using Base::Base;
};
// Primitive:
template <int kDim>
class Primitive : public Tuple<kDim> {
 public:
  // Types:
  using Base = Tuple<kDim>;
  // Constructors:
  using Base::Base;
  explicit Primitive(Base const& tuple) : Base(tuple) {}
  // Accessors and Mutators:
  Scalar const& rho() const { return this->mass; }
  Scalar const& p() const { return this->energy; }
  Scalar const& u() const { return this->momentum(0); }
  Scalar const& v() const { return this->momentum(1); }
  Scalar const& w() const { return this->momentum(2); }
  Scalar& rho() { return this->mass; }
  Scalar& p() { return this->energy; }
  Scalar& u() { return this->momentum(0); }
  Scalar& v() { return this->momentum(1); }
  Scalar& w() { return this->momentum(2); }
};
// Conservative:
template <int kDim>
struct Conservative : Tuple<kDim>{
  // Types:
  using Base = Tuple<kDim>;
  using Vector = typename Base::Vector;
  // Constructors:
  using Base::Base;
  explicit Conservative(Base const& tuple) : Base(tuple) {}
};
class IdealGas {
 private:
  static constexpr Scalar gamma_ = 1.4;

 public:
  // Constants:
  static constexpr Scalar Gamma() { return gamma_; }
  static constexpr Scalar OneOverGamma() {
    return 1 / Gamma();
  }
  static constexpr Scalar GammaPlusOne() {
    return Gamma() + 1;
  }
  static constexpr Scalar GammaPlusOneOverTwo() {
    return GammaPlusOne() / 2;
  }
  static constexpr Scalar GammaPlusOneOverFour() {
    return GammaPlusOne() / 4;
  }
  static constexpr Scalar GammaMinusOne() {
    return Gamma() - 1;
  }
  static constexpr Scalar OneOverGammaMinusOne() {
    return 1 / GammaMinusOne();
  }
  static constexpr Scalar GammaOverGammaMinusOne() {
    return Gamma() / GammaMinusOne();
  }
  static constexpr Scalar GammaMinusOneOverTwo() {
    return GammaMinusOne() / 2;
  }
  static constexpr Scalar GammaMinusOneUnderTwo() {
    return 2 / GammaMinusOne();
  }
  // Converters:
  template <int kDim>
  static Scalar GetSpeedOfSound(Primitive<kDim> const& state) {
    return state.rho() == 0 ? 0 : std::sqrt(Gamma() * state.p() / state.rho());
  }
  template <int kDim>
  static Primitive<kDim>& ConservativeToPrimitive(Tuple<kDim>* state) {
    auto& rho = state->mass;
    if (rho > 0) {
      // momentum = rho * u
      state->momentum /= rho;
      auto& u = state->momentum;
      // energy = p/(gamma - 1) + 0.5*rho*|u|^2
      state->energy -= 0.5 * rho * u.dot(u);
      state->energy *= GammaMinusOne();
      if (state->energy < 0) {
        assert(-0.0001 < state->energy);
        state->energy = 0;
      }
    } else {
      assert(rho == 0);
      state->momentum *= 0.0;
      state->energy = 0.0;
    }
    return reinterpret_cast<Primitive<kDim>&>(*state);
  }
  template <int kDim>
  static Primitive<kDim> ConservativeToPrimitive(
      Conservative<kDim> const& conservative) {
    auto primitive = Primitive{conservative};
    ConservativeToPrimitive(&primitive);
    return primitive;
  }
  template <int kDim>
  static Conservative<kDim>& PrimitiveToConservative(Tuple<kDim>* state) {
    auto& rho = state->mass;
    auto& u = state->momentum;
    // energy = p/(gamma - 1) + 0.5*rho*|u|^2
    state->energy *= OneOverGammaMinusOne();  // p / (gamma - 1)
    state->energy += 0.5 * rho * u.dot(u);  // + 0.5 * rho * |u|^2
    // momentum = rho * u
    state->momentum *= rho;
    return reinterpret_cast<Conservative<kDim>&>(*state);
  }
  template <int kDim>
  static Conservative<kDim> PrimitiveToConservative(
      Primitive<kDim> const& primitive) {
    auto conservative = Conservative{primitive};
    PrimitiveToConservative(&conservative);
    return conservative;
  }
};

}  // namespace riemann
}  // namespace buaa

#endif  //  INCLUDE_BUAA_RIEMANN_TYPES_HPP_
