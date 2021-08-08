// Copyright 2021 Weicheng Pei and Minghao Yang
#include "gtest/gtest.h"

#include "buaa/riemann/linear.hpp"

namespace buaa {
namespace riemann {

class TestLinearWaveTest : public ::testing::Test {
};
TEST_F(TestLinearWaveTest, TestFlux) {
  Scalar u_l{2.0}, u_r{1.0};
  // right running wave
  Scalar a{1.0};
  EXPECT_EQ(LinearFlux(u_l, u_r, a), LinearFlux(u_l, a));
  // left running wave
  a = -1;
  EXPECT_EQ(LinearFlux(u_l, u_r, a), LinearFlux(u_r, a));
  // standing wave
  a = 0;
  EXPECT_EQ(LinearFlux(u_l, u_r, a), LinearFlux(u_l, a));
}

}  // namespace riemann
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
