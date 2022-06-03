

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>

#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/env/geometry.hpp"
#include "libatom/utils.hpp"

using namespace otf;

std::mt19937 gen(33);
std::uniform_real_distribution<floating> dis(-10, 10);

auto main(int, char**) -> int {
  // Set up P and Q <- P

  env::Geometry<Position, Colour, Index> P;

  P.emplace_back({0, 0, 0}, 0, 0);

  for (std::size_t i = 0; i < 64; i++) {
    P.emplace_back({dis(gen), dis(gen), dis(gen)}, 0, i);
  }

  {
    auto COM = com(P);

    for (auto&& elem : P) {
      elem(Position{}) -= COM;
    }
  }

  env::Geometry Q = P;

  {  // Compute min sep
    floating delta = 9999;

    for (auto const& a : P) {
      for (auto const& b : P) {
        if (&a != &b) {
          delta = std::min(delta, norm(a(Position{}) - b(Position{})));
        }
      }
    }

    fmt::print("rmin = {}\n", delta);
  }

  // Randomly perturb P

  for (auto& elem : P) {
    elem(Position{}) += 0.01 * Vec3<floating>{dis(gen), dis(gen), dis(gen)};
  }

  {
    auto COM = com(P);

    for (auto&& elem : P) {
      elem(Position{}) -= COM;
    }
  }

  // Shuffle P

  std::shuffle(std::next(P.begin()), P.end(), gen);

  // Rotate P

  std::size_t N = 6;

  Mat3<floating> Rot{
      {std::cos(2 * M_PI / N), -std::sin(2 * M_PI / N), 0},
      {std::sin(2 * M_PI / N), +std::cos(2 * M_PI / N), 0},
      {0, 0, 1},
  };

  for (auto& elem : P) {
    elem(Position{}) = Rot * elem(Position{}).matrix();
  }

  //  Permute and transform P

  std::optional tmp = P.permute_onto(Q, 1.36);

  ASSERT(tmp, "No rotor");

  auto [rmsd, O] = *tmp;

  fmt::print("rmsd={}\n", rmsd);

  fmt::print("rmsd={}\n", env::rmsd(O, P, Q));

  for (auto&& elem : P) {
    elem(Position{}) = O * elem(Position{}).matrix();
  }

  return 0;
}