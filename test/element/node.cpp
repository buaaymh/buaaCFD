// Copyright 2021 Weicheng Pei and Minghao Yang
#include "buaa/element/node.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace element {

class NodeTest : public ::testing::Test {
 protected:
  using N2 = Node<2>;
  using N3 = Node<3>;
  const Id id{8};
  const Scalar x{1.0}, y{2.0}, z{3.0};
};
TEST_F(NodeTest, TemplateConstructor) {
  // Test N2(Real x, Real y):
  auto p2 = N2(id, x, y);
  EXPECT_EQ(p2.I(), id);
  EXPECT_EQ(p2.X(), x);
  EXPECT_EQ(p2.Y(), y);
  // Test N3(Real x, Real y, Real z):
  auto p3 = N3(id, x, y, z);
  EXPECT_EQ(p3.I(), id);
  EXPECT_EQ(p3.X(), x);
  EXPECT_EQ(p3.Y(), y);
  EXPECT_EQ(p3.Z(), z);
}

}  // namespace element
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
