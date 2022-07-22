

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <random>

#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/env/geometry.hpp"
#include "libatom/utils.hpp"

using namespace otf;

std::random_device rd{};
std::mt19937 gen(rd());
std::uniform_real_distribution<floating> dis(-1, 1);

auto main(int, char**) -> int {
  //
  env::Geometry<Position, Colour, Index> P;

  constexpr size_t N = 6;

  //   test

  {  // Set up P

    P.emplace_back({0, 0, 0}, 0, 0);

    Mat3<floating> Rot{
        {std::cos(2 * M_PI / N), -std::sin(2 * M_PI / N), 0},
        {std::sin(2 * M_PI / N), +std::cos(2 * M_PI / N), 0},
        {0, 0, 1},
    };

    Vec3<floating> x{1.0, 0.0, 0.0};

    for (size_t i = 0; i < N; i++) {
      P.emplace_back(x, 0, i + 1);
      x = Rot * x.matrix();
    }

    // Origin is the position if atom[0]
  }

  {  // Randomly perturb P

    for (size_t i = 1; i <= N; i++) {
      P[i](Position{}) += 0.05 * Vec3<floating>{dis(gen), dis(gen), dis(gen)};
    }
  }

  floating delta = 9999;

  {  // Compute min sep
    for (auto const& a : P) {
      for (auto const& b : P) {
        if (&a != &b) {
          delta = std::min(delta, norm(a(Position{}) - b(Position{})));
        }
      }
    }

    fmt::print("R_min = {:e}\n", delta);
  }

  ///

  auto primed = P;

  int count = 0;

  for_equiv_perms(P, 0.4 * delta)(1, primed, [&](Mat3<floating> const&, floating rmsd) {
    //
    fmt::print("rmsd = {:e}\n", rmsd);

    ++count;

    return false;
  });

  fmt::print("count = {}\n", count);

  return 0;
}