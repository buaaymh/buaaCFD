// Copyright 2021 Minghao Yang
#include "buaa/element/gauss.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace element {

class GaussEdgeTest : public ::testing::Test {
 protected:
  using Gauss_1 = Gauss<1>;
  using Gauss_2 = Gauss<2>;
  Gauss_1 gauss_1 = Gauss_1();
  Gauss_2 gauss_2 = Gauss_2();
};
TEST_F(GaussEdgeTest, OnePointLine) {
  EXPECT_EQ(gauss_1.x_local[0], 0);
  EXPECT_EQ(gauss_1.weights[0], 2);
  EXPECT_EQ(gauss_1.CountPoint(), 1);
}
TEST_F(GaussEdgeTest, TwoPointLine) {
  EXPECT_EQ(gauss_2.x_local[0], -0.5773502691896250);
  EXPECT_EQ(gauss_2.x_local[1], +0.5773502691896250);
  EXPECT_EQ(gauss_2.weights[0], 1);
  EXPECT_EQ(gauss_2.weights[1], 1);
  EXPECT_EQ(gauss_2.CountPoint(), 2);
}

}  // namespace element
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}