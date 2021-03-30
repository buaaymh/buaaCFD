// Copyright 2021 Weicheng Pei and Minghao Yang
#include "gtest/gtest.h"

#include "buaa/riemann/linear.hpp"

namespace buaa {
namespace riemann {

class TestLinearWaveTest : public ::testing::Test {
 protected:
  using Solver = Linear;
  using State = Solver::State;
  using Flux = Solver::Flux;
};
TEST_F(TestLinearWaveTest, TestFlux) {
  State u_l{2.0}, u_r{1.0};
  // right running wave
  auto solver = Solver(/* speed = */1.0);
  EXPECT_EQ(solver.GetFluxOnTimeAxis(u_l, u_r), solver.GetFlux(u_l));
  // left running wave
  solver = Solver(-1.0);
  EXPECT_EQ(solver.GetFluxOnTimeAxis(u_l, u_r), solver.GetFlux(u_r));
  // standing wave
  solver = Solver(0.0);
  EXPECT_EQ(solver.GetFluxOnTimeAxis(u_l, u_r), solver.GetFlux(u_l));
}

}  // namespace riemann
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
