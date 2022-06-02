

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>

#include "libatom/atom.hpp"
#include "libatom/env/geometry.hpp"

using namespace otf;

std::mt19937 gen(33);
std::uniform_real_distribution<floating> dis(-10, 10);

auto main(int, char**) -> int {
  //

  env::Geometry<Position, Colour> P;

  P.emplace_back({0, 0, 0}, 0);

  for (std::size_t i = 0; i < 64; i++) {
    P.emplace_back({dis(gen), dis(gen), dis(gen)}, 0);
  }

  auto COM = P.com();

  for (auto&& elem : P) {
    elem(Position{}) -= COM;
  }

  std::cout << P.com() << std::endl;

  ///////////////////////////////

  env::Geometry Q = P;

  Mat3<floating> Rot{
      {std::cos(1.0), -std::sin(1.0), 0},
      {std::sin(1.0), +std::cos(1.0), 0},
      {0, 0, 1},
  };

  for (auto& elem : P) {
    elem(Position{}) = Rot * elem(Position{}).matrix();
  }

  return 0;
}