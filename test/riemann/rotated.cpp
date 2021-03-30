// Copyright 2021 Weicheng Pei and Minghao Yang
#include "gtest/gtest.h"

#include "buaa/riemann/types.hpp"
#include "buaa/riemann/rotated.hpp"

namespace buaa {
namespace riemann {

class TestRotatedEulerTest : public ::testing::Test {
 protected:
  using Gas = IdealGas;
  using UnrotatedSolver = Ausm<Gas, 2>;
  using Solver = RotatedEuler<UnrotatedSolver>;
  using Vector = Solver::Vector;
  using State = Solver::State;
  Solver solver;
};
TEST_F(TestRotatedEulerTest, TestVectorConverter) {
  Vector n(+0.6, 0.8), t(-0.8, 0.6), v(3.0, 4.0), v_copy(3.0, 4.0);
  solver.Rotate(n);
  solver.GlobalToNormal(&v);
  EXPECT_EQ(v(0), v_copy.dot(n));
  EXPECT_EQ(v(1), v_copy.dot(t));
  solver.NormalToGlobal(&v);
  EXPECT_DOUBLE_EQ(v(0), v_copy(0));
  EXPECT_DOUBLE_EQ(v(1), v_copy(1));
}

}  // namespace riemann
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
