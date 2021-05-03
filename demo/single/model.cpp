// Copyright 2019 Weicheng Pei and Minghao Yang

#include <cmath>
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <string>

#include "buaa/mesh/data.hpp"
#include "buaa/mesh/dim2.hpp"
#include "buaa/mesh/vtk/reader.hpp"
#include "buaa/mesh/vtk/writer.hpp"
#include "buaa/riemann/rotated.hpp"
#include "buaa/solver/rkvr.hpp"
#include "buaa/data/path.hpp"  // defines TEST_DATA_DIR

namespace buaa {
namespace solver {

using Scalar = double;

template <class Riemann>
class SingleWaveTest {
 public:
  explicit SingleWaveTest(char** argv)
      : model_name_{argv[1]},
        mesh_name_{argv[2]},
        start_{std::atof(argv[3])},
        stop_{std::atof(argv[4])},
        n_steps_{std::atoi(argv[5])},
        output_rate_{std::atoi(argv[6])} {
    duration_ = stop_ - start_;
  }

  void Run() {
    Mesh::Cell::scalar_names.at(0) = "U";
    auto u_l = State{-1.0};
    auto u_r = State{+1.0};
    auto model = Model(model_name_);
    Riemann::global_coefficient[0] = 1.0;
    Riemann::global_coefficient[1] = 0.0;
    model.ReadMesh(test_data_dir_ + mesh_name_);
    // Set Boundary Conditions:
    constexpr auto eps = 1e-5;
    model.SetBoundaryName("left", [&](Edge& edge) {
      return std::abs(edge.Center().X() + 1.0) < eps;
    });
    model.SetBoundaryName("right", [&](Edge& edge) {
      return std::abs(edge.Center().X() - 1.0) < eps;
    });
    model.SetBoundaryName("top", [&](Edge& edge) {
      return std::abs(edge.Center().Y() - 0.05) < eps;
    });
    model.SetBoundaryName("bottom", [&](Edge& edge) {
      return std::abs(edge.Center().Y() + 0.05) < eps;
    });
    model.SetPeriodicBoundary("top", "bottom");
    model.SetPeriodicBoundary("left", "right");
    // Set Initial Conditions:
    model.SetInitialState([&](Cell& cell) {
      cell.data.u_stages[0] = cell.Integrate([&](auto point){
          return std::sin(point.X() * acos(0.0) * 4);}) / cell.Measure();
    });
    model.SetTimeSteps(duration_, n_steps_, output_rate_);
    auto output_dir = std::string("result/demo/") + model_name_;
    model.SetOutputDir(output_dir + "/");
    system(("rm -rf " + output_dir).c_str());
    system(("mkdir -p " + output_dir).c_str());
    // Commit the calculation:
    model.Calculate();
  }

 protected:
  using Jacobi = typename Riemann::Jacobi;
  using State = typename Riemann::State;
  using Flux = typename Riemann::Flux;
  static constexpr int degree = 1;
  static constexpr int num_coefficients = (degree+1) * (degree+2) / 2 - 1;
  using Stages = Eigen::Matrix<Scalar, 3, 1>;
  using Coefficients = Eigen::Matrix<Scalar, num_coefficients, 1>;
  // Types:
  struct EdgeData : public mesh::Empty {
    Scalar distance;
    Flux flux;
    Riemann riemann;
  };
  struct CellData : public mesh::Data<
      2/* dims */, 1/* scalars */, 0/* vectors */> {
   public:
    Stages u_stages;
    Coefficients coefficients;
    void Write() {
      scalars[0] = u_stages[0];
    }
  };
  using Mesh = mesh::Mesh<degree, EdgeData, CellData>;
  using Cell = typename Mesh::Cell;
  using Edge = typename Mesh::Edge;
  using Model = solver::Rkvr<Mesh, Riemann>;
  // Data:
  const std::string test_data_dir_{TEST_DATA_DIR};
  const std::string model_name_;
  const std::string mesh_name_;
  double duration_;
  double start_;
  double stop_;
  int n_steps_;
  int output_rate_;
};

}  // namespace solver
}  // namespace buaa


int main(int argc, char* argv[]) {
  if (argc == 1) {
    std::cout << "usage: model ";  // argv[0] == "model"
    std::cout << "<linear> ";
    std::cout << "<mesh> ";
    std::cout << "<start> <stop> <steps> ";
    std::cout << "<output_rate> ";
    std::cout << std::endl;
  } else if (argc == 7) {
    using Linear = buaa::riemann::RotatedLinear;
    using LinearTest = buaa::solver::SingleWaveTest<Linear>;
    if (std::strcmp(argv[1], "linear") == 0) {
      auto model = LinearTest(argv);
      model.Run();
    } else {
    }
  }
}
