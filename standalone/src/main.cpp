

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>

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

  constexpr int N = 4;

  Mat3<floating> Rot{
      {std::cos(2 * M_PI / N), -std::sin(2 * M_PI / N), 0},
      {std::sin(2 * M_PI / N), +std::cos(2 * M_PI / N), 0},
      {0, 0, 1},
  };

  Vec3<floating> basis{1, 0, 0};

  for (size_t i = 0; i < N; i++) {
    P.emplace_back(basis, 0, i + 1);
    basis = Rot * basis.matrix();
  }

  // for (std::size_t i = 0; i < 64; i++) {
  //   P.emplace_back({dis(gen), dis(gen), dis(gen)}, 0);
  // }

  {
    auto COM = P.com();

    for (auto&& elem : P) {
      elem(Position{}) -= COM;
    }
  }

  // Randomly perturb P

  // for (auto& elem : P) {
  //   elem(Position{}) += 0.005 * Vec3<floating>{dis(gen), dis(gen), dis(gen)};
  // }

  // {
  //   auto COM = P.com();

  //   for (auto&& elem : P) {
  //     elem(Position{}) -= COM;
  //   }
  // }

  for (size_t i = 0; i < 5; i++) {
    std::cout << P[i](Position{}).transpose() << std::endl;
  }

  env::Geometry Q = P;

  // Shuffle P

  // std::shuffle(std::next(P.begin()), P.end(), gen);

  // Rotate P

  // for (auto& elem : P) {
  //   elem(Position{}) = Rot * elem(Position{}).matrix();
  // }

  //  Permute and transform P

  auto [rmsd, O] = *P.permute_onto(Q, 0.15);

  for (auto&& elem : P) {
    elem(Position{}) = O * elem(Position{}).matrix();
  }

  return 0;
}