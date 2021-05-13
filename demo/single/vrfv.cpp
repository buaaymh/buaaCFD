// Copyright 2019 Weicheng Pei and Minghao Yang

#include <cmath>
#include <cstring>
#include <initializer_list>
#include <iostream>
#include <string>

#include "gtest/gtest.h"

#include "buaa/mesh/data.hpp"
#include "buaa/mesh/dim2.hpp"
#include "buaa/riemann/rotated.hpp"
#include "buaa/solver/rkvr.hpp"
#include "buaa/data/path.hpp"  // defines TEST_DATA_DIR

namespace buaa {
namespace solver {

class VrfvTest : public ::testing::Test {
 protected:
  static constexpr int degree = 1;
  static constexpr int num_coefficients = (degree+1) * (degree+2) / 2 - 1;
  // Types:
  using Stages = Eigen::Matrix<Scalar, 3, 1>;
  using Coefficients = Eigen::Matrix<Scalar, num_coefficients, 1>;
  using Riemann = buaa::riemann::RotatedLinear;
  using Jacobi = typename Riemann::Jacobi;
  using State = typename Riemann::State;
  using Flux = typename Riemann::Flux;
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
  const std::string model_name_{"vrfv"};
  const std::string mesh_name_{"simple.vtk"};
};
TEST_F(VrfvTest, Constructor) {
  Mesh::Cell::scalar_names.at(0) = "U";
  auto output_dir = std::string("result/demo/") + model_name_;
  auto model = Model(model_name_);
  model.ReadMesh(test_data_dir_ + mesh_name_);
  // Set Boundary Conditions:
  constexpr auto eps = 1e-5;
  model.SetBoundaryName("left", [&](Edge& edge) {
    return std::abs(edge.Center().X() + 0.0) < eps;
  });
  model.SetBoundaryName("right", [&](Edge& edge) {
    return std::abs(edge.Center().X() - 1.0) < eps;
  });
  model.SetBoundaryName("top", [&](Edge& edge) {
    return std::abs(edge.Center().Y() - 1.0) < eps;
  });
  model.SetBoundaryName("bottom", [&](Edge& edge) {
    return std::abs(edge.Center().Y() + 0.0) < eps;
  });
  model.SetPeriodicBoundary("top", "bottom");
  model.SetPeriodicBoundary("left", "right");
  // Set Initial Conditions:
  model.SetInitialState([&](Cell& cell) {
    cell.data.u_stages[0] = cell.Integrate([&](auto point){
        return -point.X() + point.Y();}) / cell.Measure();
  });
  model.InitializeVrMatrix();
  model.UpdataCoefficients(0);
  system(("rm -rf " + output_dir).c_str());
  system(("mkdir -p " + output_dir).c_str());
  auto filename = output_dir + std::string("/vrfv.vtu");
  auto pass = model.WriteCurrentFrame(filename);
}

}  // namespace solver
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}