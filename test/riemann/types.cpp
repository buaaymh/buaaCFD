// Copyright 2019 Weicheng Pei and Minghao Yang
#include "gtest/gtest.h"

#include "buaa/riemann/types.hpp"

namespace buaa {
namespace riemann {

class TestIdealGas : public ::testing::Test {
 protected:
  using Gas = IdealGas;
  static constexpr Scalar gamma = Gas::Gamma();
};
TEST_F(TestIdealGas, TestConverters) {
  Scalar rho{0.1}, u{+0.2}, v{-0.2}, p{0.3};
  auto primitive = Primitive<2>{rho, u, v, p};
  auto conservative = Conservative<2>{
    rho, rho*u, rho*v, p/(gamma-1) + 0.5*rho*(u*u + v*v)
  };
  EXPECT_EQ(Gas::PrimitiveToConservative(&primitive), conservative);
  auto primitive_copy = Gas::ConservativeToPrimitive(&conservative);
  EXPECT_EQ(primitive_copy.rho(), rho);
  EXPECT_DOUBLE_EQ(primitive_copy.u(), u);
  EXPECT_DOUBLE_EQ(primitive_copy.v(), v);
  EXPECT_DOUBLE_EQ(primitive_copy.p(), p);
  Primitive<2> primitive_2 = primitive_copy;
  primitive_2 += primitive_copy * 2.0;
  EXPECT_EQ(primitive_2.rho(), rho * 3);
  EXPECT_DOUBLE_EQ(primitive_2.u(), u * 3);
  EXPECT_DOUBLE_EQ(primitive_2.v(), v * 3);
  EXPECT_DOUBLE_EQ(primitive_2.p(), p * 3);
}

}  // namespace riemann
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
