#include "libatom/env/geometry.hpp"

#include <doctest/doctest.h>

#include <array>
#include <cmath>
#include <random>

#include "libatom/atom.hpp"
#include "libatom/utils.hpp"

using namespace otf;

TEST_CASE("env::geo ortho_onto, perm_onto") {
  env::Geometry<Position, Colour> P;

  P.emplace_back({0, 0, 0}, 0);

  P.emplace_back({1, 0, 0}, 0);
  P.emplace_back({0, 1, 0}, 0);
  P.emplace_back({0, 0, 1}, 0);

  env::Geometry Q = P;

  Mat3<floating> Rot{
      {std::cos(1.0), -std::sin(1.0), 0},
      {std::sin(1.0), +std::cos(1.0), 0},
      {0, 0, 1},
  };

  for (auto& elem : P) {
    elem(Position{}) = Rot * elem(Position{}).matrix();
  }

  Mat3<floating> Rp = P.ortho_onto(Q);

  floating diff = (Rp.transpose() - Rot).array().abs().sum();

  REQUIRE(diff < 0.0001);

  using std::swap;

  swap(P[1], P[2]);
  swap(P[2], P[3]);

  std::optional res = P.permute_onto(Q, 0.1);

  REQUIRE(res);

  REQUIRE((res->O.transpose() - Rot).array().abs().sum() < 0.0001);
}