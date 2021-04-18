// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_
#define INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_

#include <array>
#include <cmath>
#include <initializer_list>

#include <Eigen/Dense>

#include "buaa/element/edge.hpp"
#include "buaa/element/node.hpp"
#include "buaa/element/data.hpp"

namespace buaa {
namespace element {

class Edge;

class Triangle {
 public:
  // Types:
  using EdgeType = Edge;
  using NodeType = Node<2>;
  using PointType = Point<2>;
  using Matrix3 = Matrix<Scalar, 3, 3>;
  using Vector3 = Matrix<Scalar, 3, 1>;
  // Constructors:
  Triangle() = default;
  Triangle(Id id, const NodeType& a, const NodeType& b, const NodeType& c,
      std::initializer_list<EdgeType*> edges) : id_(id), a_(a), b_(b), c_(c), 
      edges_{*edges.begin(), *(edges.begin()+1), *(edges.begin()+2)} { 
      measure_ = GetMeasure(a, b, c);
      center_ = PointType((a + b + c) / 3);
      mat_ << 1.0, 1.0, 1.0, a.X(), b.X(), c.X(), a.Y(), b.Y(), c.Y();
  }
  // Accessors:
  Id I() const { return id_; }
  int CountVertices() const { return 3; }
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
  // Integrator:
  PointType GetGlobalXY(Scalar a, Scalar b, Scalar c) const {
    Vector3 abc = Vector3{a, b, c};
    Vector3 xy = mat_ * (abc);
    return PointType(xy(1), xy(2));
  }
  template <class Integrand>
  Scalar Integrate(Integrand&& integrand) const {
    Scalar result = 0.0;
    for (int i = 0; i < 4; ++i) {
      auto point = GetGlobalXY(local_a[i], local_b[i], local_c[i]);
      result += integrand(point.X(), point.Y()) * weights[i];
    }
    return result *= Measure();
  }
  // Geometric methods:
  Scalar Measure() const { return measure_; }
  const PointType& Center() const { return center_; }

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
  std::array<Edge*, 3> edges_;
  Scalar measure_;
  PointType center_;
  Matrix3 mat_;
};

}  // namespace element
}  // namespace buaa

#endif  //  INCLUDE_BUAA_ELEMENT_TRIANGLE_HPP_