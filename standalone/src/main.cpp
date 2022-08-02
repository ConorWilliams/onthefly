#include <fmt/core.h>
#include <fmt/ostream.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>
#include <new>
#include <optional>
#include <random>

#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/data/SoA.hpp"
#include "libatom/data/combi.hpp"
#include "libatom/data/viewSoA.hpp"
#include "libatom/env/geometry.hpp"
#include "libatom/utils.hpp"

using namespace otf;

std::random_device rd{};
std::mt19937 gen(rd());
std::uniform_real_distribution<floating> dis(-1, 1);

struct Pos : data2::MemTag<double, 3> {};

struct Cx : data2::MemTag<int> {};

struct Inertia : data2::MemTag<float, 3, 3> {};

void foo(Pos::array_cref_t) {}

auto main(int, char**) -> int {

  data2::Adaptor<Pos> p(2);

  p(Pos{}, 0) = Pos::matrix_t{1, 2, 3};

  std::cout << "p " << p[Pos{}].transpose() << std::endl;

  data2::Adaptor<Cx> q(2);

  q(Cx{}, 0);

  q(Cx{}, 0) = Cx::matrix_t{1};

  std::cout << "q " << q[Cx{}].transpose() << std::endl;

  data2::Adaptor<Pos&> v = p;
  data2::Adaptor<Pos const&> cv = p;

  std::cout << "p " << p[Pos{}].transpose() << std::endl;
  std::cout << "v " << v[Pos{}].transpose() << std::endl;
  std::cout << "cv " << cv[Pos{}].transpose() << std::endl;

  std::cout << "modify via v\n";

  v[Pos{}][0] = 99;

  std::cout << "p " << p[Pos{}].transpose() << std::endl;
  std::cout << "v " << v[Pos{}].transpose() << std::endl;
  std::cout << "cv " << cv[Pos{}].transpose() << std::endl;

  //   [[maybe_unused]] data2::Adaptor<Cx&> _ = cv;

  data2::SoA<Pos, Cx> const test(3);

  //   test[Cx{}] = 9;

  data2::SoA<Pos, Cx const&> soa{test};

  std::cout << "soa p " << soa[Pos{}].transpose() << std::endl;
  std::cout << "soa c " << soa[Cx{}].transpose() << std::endl;

  return 0;
  //

  constexpr size_t N = 6;

  //   fmt::print("align={}\n", EIGEN_MAX_ALIGN_BYTES);

  //   data::SoA<Pos, Inertia> test(2);

  //   fmt::print("size={}\n", test.size());

  //   test(Pos{}, 0) = Pos::matrix_t{1, 2, 3};
  //   test(Pos{}, 1) = Pos::matrix_t{9, 9, 9};

  //   //   SoA<Pos>

  //   //   data::SoA<Pos> t2 = test;

  //   std::cout << test[Pos{}].transpose() << std::endl;

  //   data::SoA<Pos> slice{std::move(test)};

  //   std::cout << "slice " << slice[Pos{}].transpose() << std::endl;

  //   data::ViewSoA<Pos const> view = slice;

  //   //   int i = view[Pos{}];

  //   fmt::print("modify\n");

  //   slice(Pos{}, 0)[0] = {0};
  //   //   view[Pos{}][1] = 1000;

  //   std::cout << "view  " << view[Pos{}].transpose() << std::endl;
  //   std::cout << "slice " << slice[Pos{}].transpose() << std::endl;

  return 0;

  env::Geometry<Position, Colour, Index> P;

  //   test

  { // Set up P

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

  { // Randomly perturb P

    for (size_t i = 1; i <= N; i++) {
      P[i](Position{}) += 0.05 * Vec3<floating>{dis(gen), dis(gen), dis(gen)};
    }
  }

  floating delta = 9999;

  { // Compute min sep
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