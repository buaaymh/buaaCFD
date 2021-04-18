// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_MESH_BASIS_HPP_
#define INCLUDE_BUAA_MESH_BASIS_HPP_

#include <algorithm>
#include <cmath>

#include "buaa/element/triangle.hpp"

namespace buaa {
namespace mesh {

using Scalar = element::Scalar;

template <int Order>
class Basis;

template <>
class Basis<1> {
 public:
  // Types:
  using Cell = element::Triangle;
  using Point = Cell::PointType;
  Basis(Cell* cell) : cell_{cell} {
    dx_inv_ = GetDelta(cell_->A().X(), cell_->B().X(), cell_->C().X());
    dy_inv_ = GetDelta(cell_->A().Y(), cell_->B().Y(), cell_->C().Y());
  }
  int Degree() const { return 1; }
  const Point& Center() const { return cell_->Center(); }
  Scalar Measure() const { return cell_->Measure(); }
  Scalar DxInv() const { return dx_inv_; }
  Scalar DyInv() const { return dy_inv_; }

  // Basis Functions:
  Scalar F_0_0_0(Scalar x, Scalar y) { return (x - Center().X()) * DxInv(); }
  Scalar F_0_1_0(Scalar x, Scalar y) { return DxInv(); }
  Scalar F_1_0_0(Scalar x, Scalar y) { return (y - Center().Y()) * DyInv(); }
  Scalar F_1_0_1(Scalar x, Scalar y) { return DyInv(); }

 private:
  Scalar GetDelta(Scalar a, Scalar b, Scalar c) {
    auto d = (std::max(std::max(a, b), c) - std::min(std::min(a, b), c)) * 0.5;
    return 1 / d;
  }
 private:
  Cell* cell_;
  Scalar dx_inv_;
  Scalar dy_inv_;
};

template <>
class Basis<2> : public Basis<1> {
 public:
  // Types:
  using Cell = element::Triangle;
  Basis(Cell* cell) : Basis<1>(cell) {
    auto measure_inv = 1 / Measure();
    xx_ = cell->Integrate([&](auto x, auto y){
      return std::pow(F_0_0_0(x, y), 2); }) * measure_inv;
    xy_ = cell->Integrate([&](auto x, auto y){
      return F_0_0_0(x, y) * F_1_0_0(x, y); }) * measure_inv;
    yy_ = cell->Integrate([&](auto x, auto y){
      return std::pow(F_1_0_0(x, y), 2); }) * measure_inv;
  }
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
class Basis<3> : public Basis<2> {
 public:
  // Types:
  using Cell = element::Triangle;
  Basis(Cell* cell) : Basis<2>(cell) {
    auto measure_inv = 1 / cell->Measure();
    xxx_ = cell->Integrate([&](auto x, auto y){
      return std::pow(F_0_0_0(x, y), 3); }) * measure_inv;
    xxy_ = cell->Integrate([&](auto x, auto y){
      return std::pow(F_0_0_0(x, y), 2) * F_1_0_0(x, y); }) * measure_inv;
    xyy_ = cell->Integrate([&](auto x, auto y){
      return F_0_0_0(x, y) * std::pow(F_1_0_0(x, y), 2); }) * measure_inv;
    yyy_ = cell->Integrate([&](auto x, auto y){
      return std::pow(F_1_0_0(x, y), 3); }) * measure_inv;
  }
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

}  // namespace mesh
}  // namespace buaa

#endif  //  INCLUDE_BUAA_MESH_BASIS_HPP_