// Copyright 2019 Weicheng Pei and Minghao Yang
#ifndef BUAA_SOLVER_BOUNDARY_HPP_
#define BUAA_SOLVER_BOUNDARY_HPP_

#include "buaa/mesh/dim2.hpp"

namespace buaa {
namespace solver {

using Scalar = mesh::Scalar;

template <class Mesh>
class Manager {
 public:
  // Types:
  using CellType = typename Mesh::Cell;
  using EdgeType = typename Mesh::Edge;
  using NodeType = typename Mesh::Node;
  using PointType = typename Mesh::Point;
  using Part = std::vector<EdgeType*>;
  // Construction:
  Manager() = default;
  // Mutators:
  void AddInteriorEdge(EdgeType* edge) {
    interior_edges_.emplace_back(edge);
  }
  void AddBoundaryEdge(EdgeType* edge) {
    boundary_edges_.emplace_back(edge);
  }
  template <class Visitor>
  void SetBoundaryName(std::string const& name, Visitor&& visitor) {
    name_to_part_.emplace(name, std::make_unique<Part>());
    auto& part = name_to_part_[name];
    for (auto& edge : boundary_edges_) {
      if (visitor(*edge)) {
        part->emplace_back(edge);
      }
    }
  }
  void SetPeriodicBoundary(std::string const& head, std::string const& tail) {
    periodic_part_pairs_.emplace_back(name_to_part_[head].get(),
                                      name_to_part_[tail].get());
    SetPeriodicBoundary(name_to_part_[head].get(), name_to_part_[tail].get());
  }
  void ClearBoundaryCondition() {
    if (CheckBoundaryConditions()) {
      boundary_edges_.clear();
    } else {
      throw std::length_error("Some `EdgeType`s do not have BC info.");
    }
  }
  // Iterators:
  template<class Visitor>
  void ForEachInteriorEdge(Visitor&& visit) {
    for (auto& edge : interior_edges_) {
      visit(*edge);
    }
  }
  template<class Visitor>
  void ForEachPeriodicEdge(Visitor&& visit) {
    for (auto& [left, right] : periodic_part_pairs_) {
      for (int i = 0; i < left->size(); i++) {
        visit(*(left->at(i)), *(right->at(i)));
      }
    }
  }
  void UpdatePeriodCells() {
    for (auto& [ghost, parent] : ghost_to_parent_) {
      ghost->data = parent->data;
    }
  }

 private:
  // Data members:
  std::vector<EdgeType*> interior_edges_;
  std::vector<EdgeType*> boundary_edges_;
  std::vector<std::pair<Part*, Part*>> periodic_part_pairs_;
  std::unordered_map<std::string, std::unique_ptr<Part>> name_to_part_;
  std::vector<std::pair<std::unique_ptr<CellType>, CellType*>> ghost_to_parent_;

  // Implement details:
  void SetPeriodicBoundary(Part* head, Part* tail) {
    assert(head->size() == tail->size());
    auto cmp = [](EdgeType* a, EdgeType* b) {
      auto point_a = a->Center();
      auto point_b = b->Center();
      if (point_a.Y() != point_b.Y()) {
        return point_a.Y() < point_b.Y();
      } else {
        return point_a.X() < point_b.X();
      }
    };
    std::sort(head->begin(), head->end(), cmp);
    std::sort(tail->begin(), tail->end(), cmp);
    for (int i = 0; i < head->size(); i++) {
      SewMatchingEdges(head->at(i), tail->at(i));
    }
  }
  void SewMatchingEdges(EdgeType* a, EdgeType* b) {
    auto a_positive = a->GetPositiveSide();
    auto a_negative = a->GetNegativeSide();
    auto b_positive = b->GetPositiveSide();
    auto b_negative = b->GetNegativeSide();
    CellType* cell_a = nullptr;
    CellType* cell_b = nullptr;
    auto vec_ab = PointType(b->Center() - a->Center());
    auto vec_ba = PointType(a->Center() - b->Center());
    if (a_positive == nullptr) {
      auto cell_a_ptr = MoveCell(a_negative, vec_ab);
      cell_a = cell_a_ptr.get();
      ghost_to_parent_.emplace_back(std::move(cell_a_ptr), a_negative);
    } else {
      auto cell_a_ptr = MoveCell(a_positive, vec_ab);
      cell_a = cell_a_ptr.get();
      ghost_to_parent_.emplace_back(std::move(cell_a_ptr), a_positive);
    }
    if (b_positive == nullptr) {
      auto cell_b_ptr = MoveCell(b_negative, vec_ba);
      cell_b = cell_b_ptr.get();
      ghost_to_parent_.emplace_back(std::move(cell_b_ptr), b_negative);
    } else {
      auto cell_b_ptr = MoveCell(b_positive, vec_ba);
      cell_b = cell_b_ptr.get();
      ghost_to_parent_.emplace_back(std::move(cell_b_ptr), b_positive);
    }
    if (a_positive == nullptr) { a->SetPositiveSide(cell_b); }
    else { a->SetNegativeSide(cell_b); }
    if (b_positive == nullptr) { b->SetPositiveSide(cell_a); }
    else { b->SetNegativeSide(cell_a); }
    a->data.distance = (a->GetPositiveSide()->Center() - a->GetNegativeSide()->Center()).norm();
    b->data.distance = (b->GetPositiveSide()->Center() - b->GetNegativeSide()->Center()).norm();
  }
  std::unique_ptr<CellType> MoveCell(CellType* cell, const PointType& ab_vec) {
    auto a = NodeType(cell->A().I(), cell->A().X() + ab_vec.X(), cell->A().Y() + ab_vec.Y());
    auto b = NodeType(cell->B().I(), cell->B().X() + ab_vec.X(), cell->B().Y() + ab_vec.Y());
    auto c = NodeType(cell->C().I(), cell->C().X() + ab_vec.X(), cell->C().Y() + ab_vec.Y());
    auto ab = EdgeType(a, b);
    auto bc = EdgeType(b, c);
    auto ca = EdgeType(c, a);
    return std::make_unique<CellType>(cell->I(), a, b, c, &ab, &bc, &ca);
  }
  bool CheckBoundaryConditions() {
    int n = 0;
    for (auto& [name, part] : name_to_part_) {
      n += part->size();
    }
    return n == boundary_edges_.size();
  }
};

}  // namespace solver
}  // namespace buaa

#endif  // BUAA_SOLVER_BOUNDARY_HPP_
