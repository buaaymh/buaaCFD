// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_
#define INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_

#include <Eigen/Dense>

#include "buaa/element/edge.hpp"

namespace buaa {
namespace element {

template <int kDegree>
class Triangle;

template <>
class Triangle<0> {
 public:
  // Types:
  using NodeType = Node<2>;
  using PointType = Point<2>;
  using Matrix = Eigen::Matrix<Scalar, 3, 3>;
  using Column = Eigen::Matrix<Scalar, 3, 1>;
  // Constructors:
  Triangle() = default;
  Triangle(Id id, const NodeType& a, const NodeType& b, const NodeType& c)
      : id_(id), a_(a), b_(b), c_(c) { 
      measure_ = GetMeasure(a, b, c);
      center_ = PointType((a + b + c) / 3);
      transform_mat_ << 1.0, 1.0, 1.0, a.X(), b.X(), c.X(), a.Y(), b.Y(), c.Y();
  }
  // Accessors:
  Id I() const { return id_; }
  static constexpr int CountVertices() { return 3; }
  const NodeType& A() const { return a_; }
  const NodeType& B() const { return b_; }
  const NodeType& C() const { return c_; }
  const NodeType& GetNode(int i) const {
    switch (i)  {
    case 0:
      return A();
    case 1:
      return B();
    case 2:
      return C();
    default:
      throw std::out_of_range("A `Triangle` has 3 `Point`s.");
    }
  }
  static constexpr int Degree() { return 0; }
  // Integrator:
  PointType GetGlobalXY(Scalar a, Scalar b, Scalar c) const {
    Column abc = Column{a, b, c};
    Column xy = transform_mat_ * (abc);
    return PointType(xy(1), xy(2));
  }
  template <class Value, class Integrand>
  void Integrate(Integrand&& integrand, Value* value) const {
    for (int i = 0; i < 4; ++i) {
      auto point = GetGlobalXY(local_a[i], local_b[i], local_c[i]);
      *value += integrand(point) * weights[i];
    }
    *value *= Measure();
  }
  // Geometric methods:
  Scalar Measure() const { return measure_; }
  const PointType& Center() const { return center_; }
  // Factorial
  static Scalar Factorial(int p) {
    int fac = 1;
    for (int i = 1; i <= p; ++i) { fac *= i; }
    return Scalar(fac);
  }
 private:
  static Scalar GetMeasure(const NodeType& a, const NodeType& b, const NodeType& c) {
    auto cross = (b.X() - a.X()) * (c.Y() - a.Y()) -
                 (b.Y() - a.Y()) * (c.X() - a.X());   
    return std::abs(cross) * 0.5;
  }
 private:
  static constexpr std::array<Scalar, 4> local_a{0.3333333333333333, 0.6, 0.2, 0.2};
  static constexpr std::array<Scalar, 4> local_b{0.3333333333333333, 0.2, 0.2, 0.6};
  static constexpr std::array<Scalar, 4> local_c{0.3333333333333333, 0.2, 0.6, 0.2};
  static constexpr std::array<Scalar, 4> weights{-0.5624999999999998,
                                                  0.5208333333333332,
                                                  0.5208333333333332,
                                                  0.5208333333333332};
  Id id_;
  const NodeType& a_;
  const NodeType& b_;
  const NodeType& c_;
  Scalar measure_;
  PointType center_;
  Matrix transform_mat_;
};

template <>
class Triangle<1> : public Triangle<0> {
 public:
  // Types:
  using Vector = Eigen::Matrix<Scalar, 2, 1>;
  using BasisF = Eigen::Matrix<Scalar, 2, 2>;
  Triangle() = default;
  Triangle(Id id, const NodeType& a, const NodeType& b, const NodeType& c) :
      Triangle<0>{id, a, b, c} {
    dx_inv_ = GetDelta(A().X(), B().X(), C().X());
    dy_inv_ = GetDelta(A().Y(), B().Y(), C().Y());
  }
  // Accessors:
  static constexpr int Degree() { return 1; }
  Scalar DxInv() const { return dx_inv_; }
  Scalar DyInv() const { return dy_inv_; }
  // Basis Functions:
  Scalar F_0_0_0(Scalar x, Scalar y) const { return (x - Center().X()) * DxInv(); }
  Scalar F_0_1_0(Scalar x, Scalar y) const { return DxInv(); }
  Scalar F_1_0_0(Scalar x, Scalar y) const { return (y - Center().Y()) * DyInv(); }
  Scalar F_1_0_1(Scalar x, Scalar y) const { return DyInv(); }
  // Normal Functions:
  Scalar F_0_N_1(Scalar x, Scalar y, const Scalar* n) const { return F_0_1_0(x, y) * n[0]; }
  Scalar F_1_N_1(Scalar x, Scalar y, const Scalar* n) const { return F_1_0_1(x, y) * n[1]; }
  Vector Functions(Scalar x, Scalar y) const {
    Vector result;
    result << F_0_0_0(x, y), F_1_0_0(x, y);
    return result;
  }
  BasisF GetFuncTable(const Triangle<1>& cell, const Scalar* coord, const Scalar* normal) const {
    BasisF mat = BasisF::Zero();
    mat(0, 0) = cell.F_0_0_0(coord[0], coord[1]), mat(0, 1) = cell.F_0_N_1(coord[0], coord[1], normal);
    mat(1, 0) = cell.F_1_0_0(coord[0], coord[1]), mat(1, 1) = cell.F_1_N_1(coord[0], coord[1], normal);
    return mat;
  }
 private:
  static Scalar GetDelta(Scalar a, Scalar b, Scalar c) {
    auto d = (std::max(std::max(a, b), c) - std::min(std::min(a, b), c)) * 0.5;
    return 1 / d;
  }
 private:
  Scalar dx_inv_;
  Scalar dy_inv_;
};

template <>
class Triangle<2> : public Triangle<1> {
 public:
  // Types:
  using Vector = Eigen::Matrix<Scalar, 5, 1>;
  using BasisF = Eigen::Matrix<Scalar, 5, 3>;
  Triangle() = default;
  Triangle(Id id, const NodeType& a, const NodeType& b, const NodeType& c) :
      Triangle<1>{id, a, b, c} {
    Scalar measure_inv = 1 / Measure();
    const_2 = Eigen::Matrix<Scalar, 3, 1>::Zero();
    auto func = [&](const auto& point){
      Eigen::Matrix<Scalar, 3, 1> col;
      col << std::pow(F_0_0_0(point.X(), point.Y()), 2),
             F_0_0_0(point.X(), point.Y()) * F_1_0_0(point.X(), point.Y()),
             std::pow(F_1_0_0(point.X(), point.Y()), 2);
      return col;
    };
    Integrate(func, &const_2); const_2 *= measure_inv;
  }
  // Accessors:
  static constexpr int Degree() { return 2; }
  Scalar XX() const { return const_2(0); }
  Scalar XY() const { return const_2(1); }
  Scalar YY() const { return const_2(2); }
  // Basis Functions:
  Scalar F_2_0_0(Scalar x, Scalar y) const { return std::pow(F_0_0_0(x, y), 2) - XX(); }
  Scalar F_2_1_0(Scalar x, Scalar y) const { return F_0_0_0(x, y) * DxInv() * 2; }
  Scalar F_2_2_0(Scalar x, Scalar y) const { return std::pow(DxInv(), 2) * 2; }
  // Normal Functions:
  Scalar F_2_N_1(Scalar x, Scalar y, const Scalar* n) const { return F_2_1_0(x, y) * n[0]; }
  Scalar F_2_N_2(Scalar x, Scalar y, const Scalar* n) const { return F_2_2_0(x, y) * n[0] * n[0]; }

  Scalar F_3_0_0(Scalar x, Scalar y) const { return F_0_0_0(x, y) * F_1_0_0(x, y) - XY(); }
  Scalar F_3_1_0(Scalar x, Scalar y) const { return DxInv() * F_1_0_0(x, y); }
  Scalar F_3_0_1(Scalar x, Scalar y) const { return F_0_0_0(x, y) * DyInv(); }
  Scalar F_3_1_1(Scalar x, Scalar y) const { return DxInv() * DyInv(); }
  // Normal Functions:
  Scalar F_3_N_1(Scalar x, Scalar y, const Scalar* n) const { return F_3_1_0(x, y) * n[0] + F_3_0_1(x, y) * n[1]; }
  Scalar F_3_N_2(Scalar x, Scalar y, const Scalar* n) const { return F_3_1_1(x, y) * n[0] * n[1] * 2; }

  Scalar F_4_0_0(Scalar x, Scalar y) const { return std::pow(F_1_0_0(x, y), 2) - YY(); }
  Scalar F_4_0_1(Scalar x, Scalar y) const { return F_1_0_0(x, y) * DyInv() * 2; }
  Scalar F_4_0_2(Scalar x, Scalar y) const { return std::pow(DyInv(), 2) * 2; }
  // Normal Functions:
  Scalar F_4_N_1(Scalar x, Scalar y, const Scalar* n) const { return F_4_0_1(x, y) * n[1]; }
  Scalar F_4_N_2(Scalar x, Scalar y, const Scalar* n) const { return F_4_0_2(x, y) * n[1] * n[1]; }

  Vector Functions(Scalar x, Scalar y) const {
    Vector result;
    result << F_0_0_0(x, y), F_1_0_0(x, y), F_2_0_0(x, y), F_3_0_0(x, y), F_4_0_0(x, y);
    return result;
  }
  BasisF GetFuncTable(const Triangle<2>& cell, const Scalar* coord, const Scalar* normal) const {
    BasisF mat = BasisF::Zero();
    mat(0, 0) = cell.F_0_0_0(coord[0], coord[1]), mat(0, 1) = cell.F_0_N_1(coord[0], coord[1], normal);
    mat(1, 0) = cell.F_1_0_0(coord[0], coord[1]), mat(1, 1) = cell.F_1_N_1(coord[0], coord[1], normal);
    mat(2, 0) = cell.F_2_0_0(coord[0], coord[1]), mat(2, 1) = cell.F_2_N_1(coord[0], coord[1], normal), mat(2, 2) = cell.F_2_N_2(coord[0], coord[1], normal);
    mat(3, 0) = cell.F_3_0_0(coord[0], coord[1]), mat(3, 1) = cell.F_3_N_1(coord[0], coord[1], normal), mat(3, 2) = cell.F_3_N_2(coord[0], coord[1], normal);
    mat(4, 0) = cell.F_4_0_0(coord[0], coord[1]), mat(4, 1) = cell.F_4_N_1(coord[0], coord[1], normal), mat(4, 2) = cell.F_4_N_2(coord[0], coord[1], normal);
    return mat;
  }
 private:
  Eigen::Matrix<Scalar, 3, 1> const_2;
};

template <>
class Triangle<3> : public Triangle<2> {
 public:
  // Types:
  using Vector = Eigen::Matrix<Scalar, 9, 1>;
  using BasisF = Eigen::Matrix<Scalar, 9, 4>;
  Triangle() = default;
  Triangle(Id id, const NodeType& a, const NodeType& b, const NodeType& c) :
      Triangle<2>{id, a, b, c} {
    Scalar measure_inv = 1 / Measure();
    const_3 = Eigen::Matrix<Scalar, 4, 1>::Zero();
    auto func = [&](const auto& point){
      Eigen::Matrix<Scalar, 4, 1> col;
      col << std::pow(F_0_0_0(point.X(), point.Y()), 3),
             std::pow(F_0_0_0(point.X(), point.Y()), 2) * F_1_0_0(point.X(), point.Y()),
             F_0_0_0(point.X(), point.Y()) * std::pow(F_1_0_0(point.X(), point.Y()), 2),
             std::pow(F_1_0_0(point.X(), point.Y()), 3);
      return col;
    };
    Integrate(func, &const_3); const_3 *= measure_inv;
  }
  // Accessors:
  static constexpr int Degree() { return 3; }
  Scalar XXX() const { return const_3(0); }
  Scalar XXY() const { return const_3(1); }
  Scalar XYY() const { return const_3(2); }
  Scalar YYY() const { return const_3(3); }
  // Basis Functions:
  Scalar F_5_0_0(Scalar x, Scalar y) const { return std::pow(F_0_0_0(x, y), 3) - XXX(); }
  Scalar F_5_1_0(Scalar x, Scalar y) const { return std::pow(F_0_0_0(x, y), 2) * DxInv() * 3; }
  Scalar F_5_2_0(Scalar x, Scalar y) const { return F_0_0_0(x, y) * std::pow(DxInv(), 2) * 6; }
  Scalar F_5_3_0(Scalar x, Scalar y) const { return std::pow(DxInv(), 3) * 6; }
  // Normal Functions:
  Scalar F_5_N_1(Scalar x, Scalar y, const Scalar* n) const { return F_5_1_0(x, y) * n[0]; }
  Scalar F_5_N_2(Scalar x, Scalar y, const Scalar* n) const { return F_5_2_0(x, y) * std::pow(n[0], 2); }
  Scalar F_5_N_3(Scalar x, Scalar y, const Scalar* n) const { return F_5_3_0(x, y) * std::pow(n[0], 3); }

  Scalar F_6_0_0(Scalar x, Scalar y) const { return std::pow(F_0_0_0(x, y), 2) * F_1_0_0(x, y) - XXY(); }
  Scalar F_6_1_0(Scalar x, Scalar y) const { return F_0_0_0(x, y) * F_1_0_0(x, y) * DxInv() * 2; }
  Scalar F_6_0_1(Scalar x, Scalar y) const { return std::pow(F_0_0_0(x, y), 2) * DyInv(); }
  Scalar F_6_2_0(Scalar x, Scalar y) const { return F_1_0_0(x, y) * std::pow(DxInv(), 2) * 2; }
  Scalar F_6_1_1(Scalar x, Scalar y) const { return F_0_0_0(x, y) * DxInv() * DyInv() * 2; }
  Scalar F_6_2_1(Scalar x, Scalar y) const { return std::pow(DxInv(), 2) * DyInv() * 2; }
  // Normal Functions:
  Scalar F_6_N_1(Scalar x, Scalar y, const Scalar* n) const { return F_6_1_0(x, y) * n[0] + F_6_0_1(x, y) * n[1]; }
  Scalar F_6_N_2(Scalar x, Scalar y, const Scalar* n) const { return F_6_2_0(x, y) * std::pow(n[0], 2) + F_6_1_1(x, y) * n[0] * n[1] * 2; }
  Scalar F_6_N_3(Scalar x, Scalar y, const Scalar* n) const { return F_6_2_1(x, y) * std::pow(n[0], 2) * n[1] * 3; }

  Scalar F_7_0_0(Scalar x, Scalar y) const { return F_0_0_0(x, y) * std::pow(F_1_0_0(x, y), 2) - XYY(); }
  Scalar F_7_1_0(Scalar x, Scalar y) const { return std::pow(F_1_0_0(x, y), 2) * DxInv(); }
  Scalar F_7_0_1(Scalar x, Scalar y) const { return F_0_0_0(x, y) * F_1_0_0(x, y) * DyInv() * 2; }
  Scalar F_7_0_2(Scalar x, Scalar y) const { return F_0_0_0(x, y) * std::pow(DyInv(), 2) * 2; }
  Scalar F_7_1_1(Scalar x, Scalar y) const { return F_1_0_0(x, y) * DxInv() * DyInv() * 2; }
  Scalar F_7_1_2(Scalar x, Scalar y) const { return DxInv() * std::pow(DyInv(), 2) * 2; }
  // Normal Functions:
  Scalar F_7_N_1(Scalar x, Scalar y, const Scalar* n) const { return F_7_1_0(x, y) * n[0] + F_7_0_1(x, y) * n[1]; }
  Scalar F_7_N_2(Scalar x, Scalar y, const Scalar* n) const { return F_7_1_1(x, y) * n[0] * n[1] * 2 + F_7_0_2(x, y) * std::pow(n[1], 2); }
  Scalar F_7_N_3(Scalar x, Scalar y, const Scalar* n) const { return F_7_1_2(x, y) * n[0] * std::pow(n[1], 2) * 3; }

  Scalar F_8_0_0(Scalar x, Scalar y) const { return std::pow(F_1_0_0(x, y), 3) - YYY(); }
  Scalar F_8_0_1(Scalar x, Scalar y) const { return std::pow(F_1_0_0(x, y), 2) * DyInv() * 3; }
  Scalar F_8_0_2(Scalar x, Scalar y) const { return F_1_0_0(x, y) * std::pow(DyInv(), 2) * 6; }
  Scalar F_8_0_3(Scalar x, Scalar y) const { return std::pow(DyInv(), 3) * 6; }
  // Normal Functions:
  Scalar F_8_N_1(Scalar x, Scalar y, const Scalar* n) const { return F_8_0_1(x, y) * n[1]; }
  Scalar F_8_N_2(Scalar x, Scalar y, const Scalar* n) const { return F_8_0_2(x, y) * std::pow(n[1], 2); }
  Scalar F_8_N_3(Scalar x, Scalar y, const Scalar* n) const { return F_8_0_3(x, y) * std::pow(n[1], 3); }
  Vector Functions(Scalar x, Scalar y) const {
    Vector result;
    result << F_0_0_0(x, y), F_1_0_0(x, y), F_2_0_0(x, y), F_3_0_0(x, y), F_4_0_0(x, y),
              F_5_0_0(x, y), F_6_0_0(x, y), F_7_0_0(x, y), F_8_0_0(x, y);
    return result;
  }
  BasisF GetFuncTable(const Triangle<3>& cell, const Scalar* coord, const Scalar* normal) const {
    BasisF mat = BasisF::Zero();
    // One Degree:
    mat(0, 0) = cell.F_0_0_0(coord[0], coord[1]), mat(0, 1) = cell.F_0_N_1(coord[0], coord[1], normal);
    mat(1, 0) = cell.F_1_0_0(coord[0], coord[1]), mat(1, 1) = cell.F_1_N_1(coord[0], coord[1], normal);
    // Two Degree:
    mat(2, 0) = cell.F_2_0_0(coord[0], coord[1]), mat(2, 1) = cell.F_2_N_1(coord[0], coord[1], normal), mat(2, 2) = cell.F_2_N_2(coord[0], coord[1], normal);
    mat(3, 0) = cell.F_3_0_0(coord[0], coord[1]), mat(3, 1) = cell.F_3_N_1(coord[0], coord[1], normal), mat(3, 2) = cell.F_3_N_2(coord[0], coord[1], normal);
    mat(4, 0) = cell.F_4_0_0(coord[0], coord[1]), mat(4, 1) = cell.F_4_N_1(coord[0], coord[1], normal), mat(4, 2) = cell.F_4_N_2(coord[0], coord[1], normal);
    // Three Degree:
    mat(5, 0) = cell.F_5_0_0(coord[0], coord[1]), mat(5, 1) = cell.F_5_N_1(coord[0], coord[1], normal), mat(5, 2) = cell.F_5_N_2(coord[0], coord[1], normal), mat(5, 3) = cell.F_5_N_3(coord[0], coord[1], normal);
    mat(6, 0) = cell.F_6_0_0(coord[0], coord[1]), mat(6, 1) = cell.F_6_N_1(coord[0], coord[1], normal), mat(6, 2) = cell.F_6_N_2(coord[0], coord[1], normal), mat(6, 3) = cell.F_6_N_3(coord[0], coord[1], normal);
    mat(7, 0) = cell.F_7_0_0(coord[0], coord[1]), mat(7, 1) = cell.F_7_N_1(coord[0], coord[1], normal), mat(7, 2) = cell.F_7_N_2(coord[0], coord[1], normal), mat(7, 3) = cell.F_7_N_3(coord[0], coord[1], normal);
    mat(8, 0) = cell.F_8_0_0(coord[0], coord[1]), mat(8, 1) = cell.F_8_N_1(coord[0], coord[1], normal), mat(8, 2) = cell.F_8_N_2(coord[0], coord[1], normal), mat(8, 3) = cell.F_8_N_3(coord[0], coord[1], normal);
    return mat;
  }
 private:
  Eigen::Matrix<Scalar, 4, 1> const_3;
};

}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_
