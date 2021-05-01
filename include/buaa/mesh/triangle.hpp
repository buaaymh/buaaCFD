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
  static constexpr int num_coefficients = (kDegree+1) * (kDegree+2) / 2 - 1;
 public:
  // Types:
  using Base = element::Triangle<kDegree>;
  using EdgeType = Edge<kDegree, EdgeData, CellData>;
  using NodeType = element::Node<2>;
  using PointType = element::Point<2>;
  using Matrix = Eigen::Matrix<Scalar, num_coefficients, num_coefficients>;
  using Vector = Eigen::Matrix<Scalar, num_coefficients, 1>;
  using Data = CellData;
  // Constructors:
  Triangle() = default;
  Triangle(Id id, const NodeType& a, const NodeType& b, const NodeType& c,
      EdgeType* ab, EdgeType* bc, EdgeType* ca) : Base(id, a, b, c),
      edges_{ab, bc, ca} {}
  // Accessors:
  static constexpr int CountCoefficients() { return num_coefficients; }
  // Iterators:
  template <class Visitor>
  void ForEachEdge(Visitor&& visitor) { for(auto& e : edges_) {visitor(*e);} }
  // Data:
  Data data;
  Matrix a_matrix_inv;
  Vector b_vector;
  static std::array<std::string, CellData::CountScalars()> scalar_names;
  static std::array<std::string, CellData::CountVectors()> vector_names;

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
