// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_
#define INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_

#include <array>
#include <cmath>

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
  using Matrix3 = Matrix<Scalar, 3, 3>;
  using Vector3 = Matrix<Scalar, 3, 1>;
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
    Vector3 abc = Vector3{a, b, c};
    Vector3 xy = transform_mat_ * (abc);
    return PointType(xy(1), xy(2));
  }
  template <class Integrand>
  Scalar Integrate(Integrand&& integrand) const {
    Scalar result = 0.0;
    for (int i = 0; i < 4; ++i) {
      auto point = GetGlobalXY(local_a[i], local_b[i], local_c[i]);
      result += integrand(point) * weights[i];
    }
    return result *= Measure();
  }
  // Geometric methods:
  Scalar Measure() const { return measure_; }
  const PointType& Center() const { return center_; }
  // Factorial
  Scalar Factorial(int p) {
    int fac = 1;
    for (int i = 1; i <= p; ++i) { fac *= i; }
    return Scalar(fac);
  }
 private:
  Scalar GetMeasure(const NodeType& a, const NodeType& b, const NodeType& c) {
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
 private:
  Id id_;
  const NodeType& a_;
  const NodeType& b_;
  const NodeType& c_;
  Scalar measure_;
  PointType center_;
  Matrix3 transform_mat_;
};

template <>
class Triangle<1> : public Triangle<0> {
 public:
  // Types:
  using Matrix = Eigen::Matrix<Scalar, 2, 2>;
  using Vector = Eigen::Matrix<Scalar, 2, 1>;
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
  Scalar F_0_0_0(Scalar x, Scalar y) { return (x - Center().X()) * DxInv(); }
  Scalar F_0_1_0(Scalar x, Scalar y) { return DxInv(); }
  Scalar F_1_0_0(Scalar x, Scalar y) { return (y - Center().Y()) * DyInv(); }
  Scalar F_1_0_1(Scalar x, Scalar y) { return DyInv(); }
  Vector Functions(Scalar x, Scalar y) {
    Vector result;
    result << F_0_0_0(x, y), F_1_0_0(x, y);
    return result;
  }
  // VR Initialize:
  Matrix InitializeMatWith(Scalar x, Scalar y, Triangle<1>* that, Scalar distance) {
    Matrix mat;
    Scalar p0 = std::pow(distance, -1) / std::pow(Factorial(0), 2);
    Scalar p1 = std::pow(distance, +1) / std::pow(Factorial(1), 2);
    mat(0, 0) = p0 * F_0_0_0(x, y) * that->F_0_0_0(x, y) +
                p1 * F_0_1_0(x, y) * that->F_0_1_0(x, y);
    mat(0, 1) = p0 * F_0_0_0(x, y) * that->F_1_0_0(x, y) +
                p1 * F_0_1_0(x, y) * that->F_1_0_1(x, y);
    mat(1, 0) = p0 * F_1_0_0(x, y) * that->F_0_0_0(x, y) +
                p1 * F_1_0_1(x, y) * that->F_0_1_0(x, y);
    mat(1, 1) = p0 * F_1_0_0(x, y) * that->F_1_0_0(x, y) +
                p1 * F_1_0_1(x, y) * that->F_1_0_1(x, y);
    if (I() > that->I()) { mat.transposeInPlace(); }
    return mat;
  }
  Matrix InitializeMatWith(Scalar x, Scalar y, Triangle<1>* that, Scalar distance, const PointType& ab) {
    Matrix mat;
    Scalar p0 = std::pow(distance, -1) / std::pow(Factorial(0), 2);
    Scalar p1 = std::pow(distance, +1) / std::pow(Factorial(1), 2);
    mat(0, 0) = p0 * F_0_0_0(x, y) * that->F_0_0_0(x + ab.X(), y + ab.Y()) +
                p1 * F_0_1_0(x, y) * that->F_0_1_0(x + ab.X(), y + ab.Y());
    mat(0, 1) = p0 * F_0_0_0(x, y) * that->F_1_0_0(x + ab.X(), y + ab.Y()) +
                p1 * F_0_1_0(x, y) * that->F_1_0_1(x + ab.X(), y + ab.Y());
    mat(1, 0) = p0 * F_1_0_0(x, y) * that->F_0_0_0(x + ab.X(), y + ab.Y()) +
                p1 * F_1_0_1(x, y) * that->F_0_1_0(x + ab.X(), y + ab.Y());
    mat(1, 1) = p0 * F_1_0_0(x, y) * that->F_1_0_0(x + ab.X(), y + ab.Y()) +
                p1 * F_1_0_1(x, y) * that->F_1_0_1(x + ab.X(), y + ab.Y());
    if (I() > that->I()) { mat.transposeInPlace(); }
    return mat;
  }
  Vector InitializeVecWith(Scalar x, Scalar y, Scalar distance) {
    Vector vec;
    vec(0) = F_0_0_0(x, y);
    vec(1) = F_1_0_0(x, y);
    return vec / distance;
  }

 private:
  Scalar GetDelta(Scalar a, Scalar b, Scalar c) {
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
  using Matrix = Eigen::Matrix<Scalar, 5, 5>;
  using Vector = Eigen::Matrix<Scalar, 5, 1>;
  Triangle() = default;
  Triangle(Id id, const NodeType& a, const NodeType& b, const NodeType& c) :
      Triangle<1>{id, a, b, c} {
    auto measure_inv = 1 / Measure();
    xx_ = Integrate([&](auto point){
      return std::pow(F_0_0_0(point.X(), point.Y()), 2); }) * measure_inv;
    xy_ = Integrate([&](auto point){
      return F_0_0_0(point.X(), point.Y()) * F_1_0_0(point.X(), point.Y()); }) * measure_inv;
    yy_ = Integrate([&](auto point){
      return std::pow(F_1_0_0(point.X(), point.Y()), 2); }) * measure_inv;
  }
  // Accessors:
  static constexpr int Degree() { return 2; }
  Scalar XX() const { return xx_; }
  Scalar XY() const { return xy_; }
  Scalar YY() const { return yy_; }
  // Basis Functions:
  Scalar F_2_0_0(Scalar x, Scalar y) { return std::pow(F_0_0_0(x, y), 2) - XX(); }
  Scalar F_2_1_0(Scalar x, Scalar y) { return F_0_0_0(x, y) * DxInv() * 2; }
  Scalar F_2_2_0(Scalar x, Scalar y) { return std::pow(DxInv(), 2) * 2; }

  Scalar F_3_0_0(Scalar x, Scalar y) { return F_0_0_0(x, y) * F_1_0_0(x, y) - XY(); }
  Scalar F_3_1_0(Scalar x, Scalar y) { return DxInv() * F_1_0_0(x, y); }
  Scalar F_3_0_1(Scalar x, Scalar y) { return F_0_0_0(x, y) * DyInv(); }
  Scalar F_3_1_1(Scalar x, Scalar y) { return DxInv() * DyInv(); }

  Scalar F_4_0_0(Scalar x, Scalar y) { return std::pow(F_1_0_0(x, y), 2) - YY(); }
  Scalar F_4_0_1(Scalar x, Scalar y) { return F_1_0_0(x, y) * DyInv() * 2; }
  Scalar F_4_0_2(Scalar x, Scalar y) { return std::pow(DyInv(), 2) * 2; }
 private:
  Scalar xx_;
  Scalar xy_;
  Scalar yy_;
};

template <>
class Triangle<3> : public Triangle<2> {
 public:
  // Types:
  using Matrix = Eigen::Matrix<Scalar, 9, 9>;
  using Vector = Eigen::Matrix<Scalar, 9, 1>;
  Triangle() = default;
  Triangle(Id id, const NodeType& a, const NodeType& b, const NodeType& c) :
      Triangle<2>{id, a, b, c} {
    auto measure_inv = 1 / Measure();
    xxx_ = Integrate([&](auto point){
      return std::pow(F_0_0_0(point.X(), point.Y()), 3); }) * measure_inv;
    xxy_ = Integrate([&](auto point){
      return std::pow(F_0_0_0(point.X(), point.Y()), 2) * F_1_0_0(point.X(), point.Y()); }) * measure_inv;
    xyy_ = Integrate([&](auto point){
      return F_0_0_0(point.X(), point.Y()) * std::pow(F_1_0_0(point.X(), point.Y()), 2); }) * measure_inv;
    yyy_ = Integrate([&](auto point){
      return std::pow(F_1_0_0(point.X(), point.Y()), 3); }) * measure_inv;
  }
  // Accessors:
  static constexpr int Degree() { return 3; }
  Scalar XXX() const { return xxx_; }
  Scalar XXY() const { return xxy_; }
  Scalar XYY() const { return xyy_; }
  Scalar YYY() const { return yyy_; }
  // Basis Functions:
  Scalar F_5_0_0(Scalar x, Scalar y) { return std::pow(F_0_0_0(x, y), 3) - XXX(); }
  Scalar F_5_1_0(Scalar x, Scalar y) { return std::pow(F_0_0_0(x, y), 2) * DxInv() * 3; }
  Scalar F_5_2_0(Scalar x, Scalar y) { return F_0_0_0(x, y) * std::pow(DxInv(), 2) * 6; }
  Scalar F_5_3_0(Scalar x, Scalar y) { return std::pow(DxInv(), 3) * 6; }

  Scalar F_6_0_0(Scalar x, Scalar y) { return std::pow(F_0_0_0(x, y), 2) * F_1_0_0(x, y) - XXY(); }
  Scalar F_6_1_0(Scalar x, Scalar y) { return F_0_0_0(x, y) * F_1_0_0(x, y) * DxInv() * 2; }
  Scalar F_6_2_0(Scalar x, Scalar y) { return F_1_0_0(x, y) * std::pow(DxInv(), 2) * 2; }
  Scalar F_6_0_1(Scalar x, Scalar y) { return std::pow(F_0_0_0(x, y), 2) * DyInv(); }
  Scalar F_6_1_1(Scalar x, Scalar y) { return F_0_0_0(x, y) * DxInv() * DyInv() * 2; }
  Scalar F_6_2_1(Scalar x, Scalar y) { return std::pow(DxInv(), 2) * DyInv() * 2; }

  Scalar F_7_0_0(Scalar x, Scalar y) { return F_0_0_0(x, y) * std::pow(F_1_0_0(x, y), 2) - XYY(); }
  Scalar F_7_1_0(Scalar x, Scalar y) { return std::pow(F_1_0_0(x, y), 2) * DxInv(); }
  Scalar F_7_0_1(Scalar x, Scalar y) { return F_0_0_0(x, y) * F_1_0_0(x, y) * DyInv() * 2; }
  Scalar F_7_0_2(Scalar x, Scalar y) { return F_0_0_0(x, y) * std::pow(DyInv(), 2) * 2; }
  Scalar F_7_1_1(Scalar x, Scalar y) { return F_1_0_0(x, y) * DxInv() * DyInv() * 2; }
  Scalar F_7_1_2(Scalar x, Scalar y) { return DxInv() * std::pow(DyInv(), 2) * 2; }

  Scalar F_8_0_0(Scalar x, Scalar y) { return std::pow(F_1_0_0(x, y), 3) - YYY(); }
  Scalar F_8_0_1(Scalar x, Scalar y) { return std::pow(F_1_0_0(x, y), 2) * DyInv() * 3; }
  Scalar F_8_0_2(Scalar x, Scalar y) { return F_1_0_0(x, y) * std::pow(DyInv(), 2) * 6; }
  Scalar F_8_0_3(Scalar x, Scalar y) { return std::pow(DyInv(), 3) * 6; }
 private:
  Scalar xxx_;
  Scalar xxy_;
  Scalar xyy_;
  Scalar yyy_;
};


}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_