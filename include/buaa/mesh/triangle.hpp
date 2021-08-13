// Copyright 2021 Minghao Yang
#ifndef INCLUDE_BUAA_MESH_TRIANGLE_HPP_
#define INCLUDE_BUAA_MESH_TRIANGLE_HPP_

#include <array>
#include <string>

#include <Eigen/Dense>

#include "buaa/element/triangle.hpp"
#include "buaa/mesh/edge.hpp"
#include "buaa/mesh/data.hpp"

namespace buaa {
namespace mesh {

template <int kDegree, class EdgeData, class CellData>
class Edge;

template <int   kDegree,
          class EdgeData = Empty,
          class CellData = Empty>
class Triangle;

template <int kDegree, class EdgeData, class CellData>
class Triangle : public element::Triangle<kDegree> {
 private:
  static constexpr int nCoef = (kDegree+1) * (kDegree+2) / 2 - 1;
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // Types:
  using Base = element::Triangle<kDegree>;
  using CellType = Triangle<kDegree, EdgeData, CellData>;
  using EdgeType = Edge<kDegree, EdgeData, CellData>;
  using NodeType = element::Node<2>;
  using PointType = element::Point<2>;
  using Matrix = Eigen::Matrix<Scalar, nCoef, nCoef>;
  using Vector = Eigen::Matrix<Scalar, nCoef, 1>;
  using Matrix3V = Eigen::Matrix<Scalar, nCoef, 3>;
  using BasisF = Eigen::Matrix<Scalar, nCoef, kDegree+1>;
  using Data = CellData;
  // Constructors:
  Triangle() = default;
  Triangle(Id id, const NodeType& a, const NodeType& b, const NodeType& c,
      EdgeType* ab, EdgeType* bc, EdgeType* ca) : Base(id, a, b, c),
      edges_{ab, bc, ca} {}
  Triangle(const CellType&) = default;
  Triangle& operator=(const Triangle&);
  // Accessors:
  static constexpr int CountCoef() { return nCoef; }
  bool Contains(const EdgeType* edge) const {
    for (int i = 0; i < 3; ++i) {
      if (edges_[i] == edge) { return true; }
    }
    return false;
  }
  // Iterators:
  template <class Visitor>
  void ForEachEdge(Visitor&& visitor) { for(auto& e : edges_) {visitor(*e);} }
  // Initialize VR Matrix and Vector:
  void InitializeAmatInv() {
    Matrix a_matrix = Matrix::Zero();
    ForEachEdge([&](EdgeType& edge) {
      Matrix temp = Matrix::Zero();
      Scalar normal[2] = {edge.GetNormalX(), edge.GetNormalY()};
      edge.Integrate([&](const PointType& point) {
        return GetMatAt(point.X(), point.Y(), *this, edge.distance, normal);
      }, &temp);
      a_matrix += temp;
    });
    a_matrix_inv = a_matrix.inverse();
  }
  void InitializeBvecMat() {
    b_vector_mat = Matrix3V::Zero();
    for (int i = 0; i < 3; ++i) {
      Vector temp = Vector::Zero();
      edges_[i]->Integrate([&](const PointType& point) {
        return GetVecAt(point.X(), point.Y(), edges_[i]->distance);
      }, &temp);
      b_vector_mat.col(i) = temp;
    }
  }
  static void GetPArray(Scalar distance, int degree, Scalar* p) {
   for (int i = 0; i <= degree; ++i)
      p[i] = std::pow(distance, 2*i-1) / std::pow(Factorial(i), 2);
  }
  Matrix GetMatAt(Scalar x, Scalar y, const CellType& that, Scalar distance,
                  Scalar* n) const {
    Scalar p[kDegree+1]; GetPArray(distance, kDegree, p); Scalar coord[] = {x, y};
    // this cell
    BasisF i = GetFuncTable(*this, coord, n);
    // that cell
    BasisF j = GetFuncTable(that, coord, n);
    Matrix mat = Matrix::Zero();
    for (int m = 0; m != nCoef; ++m) {
      for (int n = 0; n != nCoef; ++n) {
        for (int k = 0; k != kDegree+1; ++k) mat(m, n) += p[k] * i(n, k) * j(m, k);
      }
    }
    if (Base::I() > that.I()) { mat.transposeInPlace(); }
    return mat;
  }
  Matrix GetMatAt(Scalar x, Scalar y, const CellType& that, Scalar distance,
                  const PointType& ab, Scalar* n) const {
    Scalar p[kDegree+1]; GetPArray(distance, kDegree, p);
    Scalar coord[] = {x, y}, coord_ab[] = {x + ab.X(), y + ab.Y()};
    // this cell
    BasisF i = GetFuncTable(*this, coord, n);
    // that cell
    BasisF j = GetFuncTable(that, coord_ab, n);
    Matrix mat = Matrix::Zero();
    for (int m = 0; m != nCoef; ++m) {
      for (int n = 0; n != nCoef; ++n) {
        for (int k = 0; k != kDegree+1; ++k) mat(m, n) += p[k] * i(n, k) * j(m, k);
      }
    }
    if (Base::I() > that.I()) { mat.transposeInPlace(); }
    return mat;
  }
  Vector GetVecAt(Scalar x, Scalar y, Scalar distance) const {
    return Functions(x, y) / distance;
  }
  // Polynomial:
  Scalar Polynomial(const PointType& point) const {
    return this->Functions(point.X(), point.Y()).transpose() * data.coefficients;
  }
  // Data:
  static std::array<std::string, CellData::CountScalars()> scalar_names;
  static std::array<std::string, CellData::CountVectors()> vector_names;
  Data data;
  Matrix a_matrix_inv;
  Vector b_vector;
  Matrix3V b_vector_mat;

 private:
  std::array<EdgeType*, 3> edges_;
};

template <int kDegree, class EdgeData, class CellData>
std::array<std::string, CellData::CountScalars()>
Triangle<kDegree, EdgeData, CellData>::scalar_names;

template <int kDegree, class EdgeData, class CellData>
std::array<std::string, CellData::CountVectors()>
Triangle<kDegree, EdgeData, CellData>::vector_names;

}  // namespace mesh
}  // namespace buaa

#endif  //  INCLUDE_BUAA_MESH_TRIANGLE_HPP_
