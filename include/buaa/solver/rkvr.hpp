// Copyright 2021 Minghao Yang
#ifndef BUAA_SOLVER_RKVR_HPP_
#define BUAA_SOLVER_RKVR_HPP_

#include <cmath>
#include <cstdio>
#include <memory>
#include <set>
#include <string>

#include "buaa/mesh/dim2.hpp"
#include "buaa/mesh/vtk/reader.hpp"
#include "buaa/mesh/vtk/writer.hpp"
#include "buaa/solver/boundary.hpp"

namespace buaa {
namespace solver {

template <class Mesh, class Riemann>
class Rkvr {
  using EdgeType = typename Mesh::Edge;
  using CellType = typename Mesh::Cell;
  using State = typename Riemann::State;
  using FluxType = typename Riemann::Flux;
  using Reader = mesh::vtk::Reader<Mesh>;
  using Writer = mesh::vtk::Writer<Mesh>;
  static constexpr int degree = Mesh::kDegree;

 public:
  explicit Rkvr(std::string const& name) : model_name_(name) {}
  bool ReadMesh(std::string const& file_name) {
    reader_ = Reader();
    if (reader_.ReadFromFile(file_name)) {
      mesh_ = reader_.GetMesh();
      Preprocess();
      return true;
    } else {
      return false;
    }
  }
  // Mutators:
  template <class Visitor>
  void SetInitialState(Visitor&& visitor) {
    mesh_->ForEachCell(visitor);
  }
  void SetTimeSteps(double duration, int n_steps, int refresh_rate) {
    duration_ = duration;
    n_steps_ = n_steps;
    step_size_ = duration / n_steps;
    refresh_rate_ = refresh_rate;
  }
  void SetOutputDir(std::string dir) {
    dir_ = dir;
  }
  template <class Visitor>
  void SetBoundaryName(std::string const& name, Visitor&& visitor) {
    edge_manager_.SetBoundaryName(name, visitor);
  }
  void SetPeriodicBoundary(std::string const& name_a,
                           std::string const& name_b) {
    edge_manager_.SetPeriodicBoundary(name_a, name_b);
  }
  // Major computation:
  void Calculate() {
    edge_manager_.ClearBoundaryCondition();
    writer_ = Writer();
    // Write the frame of initial state:
    auto filename = dir_ + model_name_ + "." + std::to_string(0) + ".vtu";
    bool pass = WriteCurrentFrame(filename);
    assert(pass);
    InitializeVrMatrix();
    for (int i = 1; i <= n_steps_ && pass; i++) {
      // Runge-Kutta three steps :
      RungeKutta3Stepper();
      if (i % refresh_rate_ == 0) {
        filename = dir_ + model_name_ + "." +std::to_string(i) + ".vtu";
        pass = WriteCurrentFrame(filename);
        std::printf("Progress: %d/%d\n", i, n_steps_);
      }
    }
  }

 private:
  bool WriteCurrentFrame(std::string const& filename) {
    mesh_->ForEachCell([&](CellType& cell) {
      cell.data.Write();
    });
    writer_.SetMesh(mesh_.get());
    return writer_.WriteToFile(filename);
  }
  void Preprocess() {
    mesh_->ForEachEdge([&](EdgeType& edge){
      auto length = edge.Measure();
      auto n1 = (edge.Tail().Y() - edge.Head().Y()) / length;
      auto n2 = (edge.Head().X() - edge.Tail().X()) / length;
      edge.data.riemann.Rotate(n1, n2);
      auto cell_l = edge.GetPositiveSide();
      auto cell_r = edge.GetNegativeSide();
      if (cell_l && cell_r) {
        edge_manager_.AddInteriorEdge(&edge);
        edge.data.distance = (cell_l->Center() - cell_r->Center()).norm();
      } else {
        edge_manager_.AddBoundaryEdge(&edge);
      }
    });
  }
  void RungeKutta3Stepper() {
    GetFluxOnEachEdge(0);
    mesh_->ForEachCell([&](CellType& cell) {
      auto rhs = GetRHS(cell);
      cell.data.u_stages[1] = cell.data.u_stages[0] +
                              rhs * step_size_ / cell.Measure();
    });
    GetFluxOnEachEdge(1);
    mesh_->ForEachCell([&](CellType& cell) {
      auto rhs = GetRHS(cell);
      cell.data.u_stages[2] = cell.data.u_stages[0] * 0.75 + 
                             (cell.data.u_stages[1] +
                              rhs * step_size_ / cell.Measure()) * 0.25;
    });
    GetFluxOnEachEdge(2);
    mesh_->ForEachCell([&](CellType& cell) {
      auto rhs = GetRHS(cell);
      cell.data.u_stages[0] = cell.data.u_stages[0] / 3 + 
                             (cell.data.u_stages[2] +
                              rhs * step_size_ / cell.Measure()) * 2 / 3;
    });
  }
  void GetFluxOnEachEdge(int stage) {
    // UpdataCoefficients(stage);
    mesh_->ForEachEdge([&](EdgeType& edge) {
      auto& riemann_ = edge.data.riemann;
      auto cell_l = edge.GetPositiveSide();
      auto cell_r = edge.GetNegativeSide();
      auto const& u_l = cell_l->data.u_stages[stage];
      auto const& u_r = cell_r->data.u_stages[stage];
      edge.data.flux = riemann_.GetFluxOnTimeAxis(u_l, u_r);
      edge.data.flux *= edge.Measure();
    });
  }

  FluxType GetRHS(CellType& cell) {
    auto rhs = FluxType();
    cell.ForEachEdge([&](EdgeType& edge) {
      if (edge.GetPositiveSide() == &cell) {
        rhs -= edge.data.flux;
      } else {
        rhs += edge.data.flux;
      }
    });
    return rhs;
  }
  void InitializeVrMatrix() {
    mesh_->ForEachCell([&](CellType& cell) {
      // cell.InitializeAinv();
    });
    mesh_->ForEachEdge([&](EdgeType& edge) {
      // edge.InitializeBmat();
    });
  }
  void UpdataCoefficients(int stage) {
    for (int i = 0; i < 5; ++i) {
      mesh_->ForEachCell([&](CellType& cell) {
        cell.data.coefficients = -0.3 * cell.data.coefficients;
        cell.ForEachEdge([&](CellType& edge) {
          auto neighber = edge.GetOpposite(&cell);
          if (neighber->I() < cell.I()) {
            cell.data.coefficients += cell.a_matrix_inv * edge.b_matrix *
                neighber->data.coefficients * 1.3;
          } else {
            cell.data.coefficients += cell.a_matrix_inv * edge.b_matrix.transpose() *
                neighber->data.coefficients * 1.3;
          }
          cell.data.coefficients += cell.a_matrix_inv * cell.b_vector * 1.3 *
              (neighber->data.u_stages[stage] - cell.data.u_stages[stage]);
        });
      });
    }
  }
  std::string model_name_;
  Reader reader_;
  Writer writer_;
  std::unique_ptr<Mesh> mesh_;
  double duration_;
  int n_steps_;
  double step_size_;
  std::string dir_;
  int refresh_rate_;
  std::set<EdgeType*> inside_edge_;
  Manager<Mesh> edge_manager_;
};

}  // namespace solver
}  // namespace buaa

#endif  // BUAA_SOLVER_RKVR_HPP_
