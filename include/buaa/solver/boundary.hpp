// Copyright 2019 Weicheng Pei and Minghao Yang
#ifndef BUAA_SOLVER_BOUNDARY_HPP_
#define BUAA_SOLVER_BOUNDARY_HPP_

#include <omp.h>

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
    #pragma omp parallel for
    for (auto& edge : interior_edges_) {
      visit(*edge);
    }
  }
  template<class Visitor>
  void ForEachPeriodicEdge(Visitor&& visit) {
    for (auto& [left, right] : periodic_part_pairs_) {
      #pragma omp parallel for
      for (int i = 0; i < left->size(); i++) {
        visit(*(left->at(i)), *(right->at(i)));
      }
    }
  }

 private:
  // Data members:
  std::vector<EdgeType*> interior_edges_;
  std::vector<EdgeType*> boundary_edges_;
  std::vector<std::pair<Part*, Part*>> periodic_part_pairs_;
  std::unordered_map<std::string, std::unique_ptr<Part>> name_to_part_;

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
    auto dist_ab = PointType(a->Center() - b->Center());
    if (a_positive == nullptr) {
      if (b_positive == nullptr) {
        a->SetPositiveSide(b_negative);
        b->SetPositiveSide(a_negative);
        dist_ab -= a_negative->Center() - b_negative->Center();
      } else {
        a->SetPositiveSide(b_positive);
        b->SetNegativeSide(a_negative);
        dist_ab -= a_negative->Center() - b_positive->Center();
      }
    } else {
      if (b_positive == nullptr) {
        a->SetNegativeSide(b_negative);
        b->SetPositiveSide(a_positive);
        dist_ab -= a_positive->Center() - b_negative->Center();
      } else {
        a->SetNegativeSide(b_positive);
        b->SetNegativeSide(a_positive);
        dist_ab -= a_positive->Center() - b_positive->Center();
      }
    }
    a->distance = dist_ab.norm();
    b->distance = dist_ab.norm();
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
