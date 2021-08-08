// Copyright 2021 Weicheng Pei and Minghao Yang
#include "gtest/gtest.h"

#include "buaa/riemann/ausm.hpp"

namespace buaa {
namespace riemann {

class Ausm2dTest : public ::testing::Test {
 protected:
  using Solver = Ausm<IdealGas, 2>;
  using State = Solver::State;
  using Flux = Solver::FluxType;
  using Vector = Solver::Vector;
  Solver solver;
  Scalar v__left{1.5}, v_right{2.5};
  static void CompareFlux(Flux const& lhs, Flux const& rhs) {
    EXPECT_EQ(lhs.mass, rhs.mass);
    EXPECT_EQ(lhs.energy, rhs.energy);
    EXPECT_EQ(lhs.momentum(0), rhs.momentum(0));
    EXPECT_EQ(lhs.momentum(1), rhs.momentum(1));
  }
};
TEST_F(Ausm2dTest, TestSod) {
  State  left{1.000, 0.0, v__left, 1.0};
  State right{0.125, 0.0, v_right, 0.1};
  CompareFlux(solver.GetFlux(left, right),
              solver.GetFlux({0.426319, +0.927453, v__left, 0.303130}));
  CompareFlux(solver.GetFlux(right, left),
              solver.GetFlux({0.426319, -0.927453, v__left, 0.303130}));
}
TEST_F(Ausm2dTest, TestShockCollision) {
  State  left{5.99924, 19.5975, v__left, 460.894};
  State right{5.99242, 6.19633, v_right, 46.0950};
  CompareFlux(solver.GetFlux(left, right),
              solver.GetFlux({5.99924, 19.5975, v__left, 460.894}));
}
TEST_F(Ausm2dTest, TestBlastFromLeft) {
  State  left{1.0, 0.0, v__left, 1e+3};
  State right{1.0, 0.0, v_right, 1e-2};
  CompareFlux(solver.GetFlux(left, right),
              solver.GetFlux({0.575062, 19.59745, v__left, 460.8938}));
}
TEST_F(Ausm2dTest, TestBlastFromRight) {
  State  left{1.0, 0.0, v__left, 1e-2};
  State right{1.0, 0.0, v_right, 1e+2};
  CompareFlux(solver.GetFlux(left, right),
              solver.GetFlux({0.575113, -6.196328, v_right, 46.09504}));
}
TEST_F(Ausm2dTest, TestAlmostVaccumed) {
  State  left{1.0, -2.0, v__left, 0.4};
  State right{1.0, +2.0, v_right, 0.4};
  CompareFlux(solver.GetFlux(left, right),
              solver.GetFlux({0.21852, 0.0, v__left, 0.001894}));
  CompareFlux(solver.GetFlux(left, right),
              solver.GetFlux({0.21852, 0.0, v_right, 0.001894}));
}
TEST_F(Ausm2dTest, TestVaccumed) {
  State  left{1.0, -4.0, v__left, 0.4};
  State right{1.0, +4.0, v_right, 0.4};
  CompareFlux(solver.GetFlux(left, right),
              solver.GetFlux({0.0, 0.0, v__left, 0.0}));
  CompareFlux(solver.GetFlux(left, right),
              solver.GetFlux({0.0, 0.0, v_right, 0.0}));
}
TEST_F(Ausm2dTest, TestNormal) {
  Vector n(+0.6, 0.8), t(-0.8, 0.6), v(3.0, 4.0), v_copy(3.0, 4.0);
  solver.GlobalToNormal(&v, n);
  EXPECT_EQ(v(0), v_copy.dot(n));
  EXPECT_EQ(v(1), v_copy.dot(t));
  solver.NormalToGlobal(&v, n);
  EXPECT_EQ(v(0), v_copy(0));
  EXPECT_EQ(v(1), v_copy(1));
}

}  // namespace riemann
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
