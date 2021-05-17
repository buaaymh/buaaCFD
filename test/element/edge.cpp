// Copyright 2021 Minghao Yang
#include "buaa/element/edge.hpp"
#include "buaa/element/triangle.hpp"

#include "gtest/gtest.h"

namespace buaa {
namespace element {

class EdgeTest : public ::testing::Test {
 protected:
  using E0 = Edge<0>;
  using E1 = Edge<1>;
  using E2 = Edge<2>;
  using E3 = Edge<3>;
  using Point = element::Point<2>;
  using Node  = element::Node<2>;
  Node head{0, 0.3, 0.0}, tail{1, 0.0, 0.4};
};
TEST_F(EdgeTest, Constructor) {
  auto e0 = E0(head, tail);
  EXPECT_EQ(e0.Head().X(), head.X());
  EXPECT_EQ(e0.Tail().X(), tail.X());
  EXPECT_EQ(e0.Head().Y(), head.Y());
  EXPECT_EQ(e0.Tail().Y(), tail.Y());
  EXPECT_EQ(e0.Measure(), 0.5);
  EXPECT_EQ(e0.Center().X() * 2, head.X() + tail.X());
  EXPECT_EQ(e0.Center().Y() * 2, head.Y() + tail.Y());
  auto e1 = E1(head, tail);
  EXPECT_EQ(e1.Head().X(), head.X());
  EXPECT_EQ(e1.Tail().X(), tail.X());
  EXPECT_EQ(e1.Head().Y(), head.Y());
  EXPECT_EQ(e1.Tail().Y(), tail.Y());
  EXPECT_EQ(e1.Measure(), 0.5);
  EXPECT_EQ(e1.Center().X() * 2, head.X() + tail.X());
  EXPECT_EQ(e1.Center().Y() * 2, head.Y() + tail.Y());
  auto e2 = E2(head, tail);
  EXPECT_EQ(e2.Head().X(), head.X());
  EXPECT_EQ(e2.Tail().X(), tail.X());
  EXPECT_EQ(e2.Head().Y(), head.Y());
  EXPECT_EQ(e2.Tail().Y(), tail.Y());
  EXPECT_EQ(e2.Measure(), 0.5);
  EXPECT_EQ(e2.Center().X() * 2, head.X() + tail.X());
  EXPECT_EQ(e2.Center().Y() * 2, head.Y() + tail.Y());
  auto e4 = E3(head, tail);
  EXPECT_EQ(e4.Head().X(), head.X());
  EXPECT_EQ(e4.Tail().X(), tail.X());
  EXPECT_EQ(e4.Head().Y(), head.Y());
  EXPECT_EQ(e4.Tail().Y(), tail.Y());
  EXPECT_EQ(e4.Measure(), 0.5);
  EXPECT_EQ(e4.Center().X() * 2, head.X() + tail.X());
  EXPECT_EQ(e4.Center().Y() * 2, head.Y() + tail.Y());
}
TEST_F(EdgeTest, Methods) {
  auto e0 = E0(head, tail);
  EXPECT_EQ(e0.Degree(), 0);
  EXPECT_EQ(e0.CountQuadPoints(), 1);
  EXPECT_EQ(e0.GetGauss().x_local[0], 0);
  EXPECT_EQ(e0.GetGauss().weights[0], 2);
  auto e1 = E1(head, tail);
  EXPECT_EQ(e1.Degree(), 1);
  EXPECT_EQ(e1.CountQuadPoints(), 1);
  EXPECT_EQ(e1.GetGauss().x_local[0], 0);
  EXPECT_EQ(e1.GetGauss().weights[0], 2);
  auto e2 = E2(head, tail);
  EXPECT_EQ(e2.Degree(), 2);
  EXPECT_EQ(e2.CountQuadPoints(), 2);
  EXPECT_EQ(e2.GetGauss().x_local[0], -0.5773502691896250);
  EXPECT_EQ(e2.GetGauss().x_local[1], +0.5773502691896250);
  EXPECT_EQ(e2.GetGauss().weights[0], 1);
  EXPECT_EQ(e2.GetGauss().weights[1], 1);
  auto e3 = E3(head, tail);
  EXPECT_EQ(e3.Degree(), 3);
  EXPECT_EQ(e3.CountQuadPoints(), 2);
  EXPECT_EQ(e3.GetGauss().x_local[0], -0.5773502691896250);
  EXPECT_EQ(e3.GetGauss().x_local[1], +0.5773502691896250);
  EXPECT_EQ(e3.GetGauss().weights[0], 1);
  EXPECT_EQ(e3.GetGauss().weights[1], 1);
}
TEST_F(EdgeTest, GaussPoints) {
  auto e0 = E0(head, tail);
  e0.ForEachQuadPoint([&](Point const& point) {
    EXPECT_EQ(point.X(), e0.Center().X());
    EXPECT_EQ(point.Y(), e0.Center().Y());
  });
  auto e1 = E1(head, tail);
  e1.ForEachQuadPoint([&](Point const& point) {
    EXPECT_EQ(point.X(), e1.Center().X());
    EXPECT_EQ(point.Y(), e1.Center().Y());
  });
  auto e2 = E2(head, tail);
  auto integrand = 0.0;
  int i = 0;
  e2.ForEachQuadPoint([&](Point const& point) {
    integrand += point.X() * e2.GetGauss().weights[i++];
  });
  EXPECT_EQ(integrand*0.5*e2.Measure(), 0.075);
  integrand = e2.Integrate([](const Point& point) { return point.X(); });
  EXPECT_EQ(integrand, 0.075);
  integrand = 0.0, i = 0;
  e2.ForEachQuadPoint([&](Point const& point) {
    integrand += point.Y() * e2.GetGauss().weights[i++];
  });
  EXPECT_EQ(integrand*0.5*e2.Measure(), 0.1);
  integrand = e2.Integrate([](const Point& point) { return point.Y(); });
  EXPECT_EQ(integrand, 0.1);
  // Three Degree
  auto eps = 10e-6;
  Node orgin{2, 0.0, 0.0};
  auto e3 = E3(orgin, head);
  integrand = e3.Integrate([](const Point& point) { return std::pow(point.X(), 3); });
  EXPECT_NEAR(integrand, 0.002025, eps);
  auto e4 = E3(orgin, tail);
  integrand = e4.Integrate([](const Point& point) { return std::pow(point.Y(), 3); });
  EXPECT_NEAR(integrand, 0.0064, eps);
  // Matrix Integrand
  using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
  auto mat = Matrix2();
  mat << 1, 0, 0, 1;
  auto mat_integrand = e1.IntegrateM([&](const Point& point) { return mat; });
  EXPECT_NEAR(mat_integrand(0,0), e1.Measure(), eps);
  EXPECT_NEAR(mat_integrand(1,0), 0.0, eps);
  EXPECT_NEAR(mat_integrand(0,1), 0.0, eps);
  EXPECT_NEAR(mat_integrand(1,1), e1.Measure(), eps);
}

}  // namespace element
}  // namespace buaa

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
