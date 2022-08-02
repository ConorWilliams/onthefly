#define EIGEN_RUNTIME_NO_MALLOC

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

template <typename T>
void pretty_print(std::string_view v, T const& x) {
  Eigen::internal::set_is_malloc_allowed(true);
  std::cout << v << ':' << x.transpose() << std::endl;
  Eigen::internal::set_is_malloc_allowed(false);
}

auto main(int, char**) -> int {

  Eigen::internal::set_is_malloc_allowed(true);

  data2::detail::Adaptor<Pos> p(2);

  Eigen::internal::set_is_malloc_allowed(false);

  fmt::print("a\n");

  p(Pos{}, 0) = Pos::matrix_t{1, 2, 3};

  fmt::print("b\n");

  p[Pos{}].transpose();

  fmt::print("b2\n");

  pretty_print("p", p[Pos{}]);

  fmt::print("b3\n");

  Eigen::internal::set_is_malloc_allowed(true);

  data2::detail::Adaptor<Cx> q(2);

  Eigen::internal::set_is_malloc_allowed(false);

  fmt::print("c\n");

  q(Cx{}, 0);

  fmt::print("d\n");

  q(Cx{}, 0) = Cx::matrix_t{1};

  pretty_print("q", q[Cx{}]);

  data2::detail::Adaptor<Pos&> v = p;
  data2::detail::Adaptor<Pos const&> cv = p;

  pretty_print("p", p[Pos{}]);
  pretty_print("v", v[Pos{}]);
  pretty_print("cv", cv[Pos{}]);

  std::cout << "modify via v\n";

  v[Pos{}][0] = 99;

  pretty_print("p", p[Pos{}]);
  pretty_print("v", v[Pos{}]);
  pretty_print("cv", cv[Pos{}]);

  Eigen::internal::set_is_malloc_allowed(true);
  data2::SoA<Pos, Cx> test(2);
  Eigen::internal::set_is_malloc_allowed(false);

  test(Cx{}, 0) = 69;

  test[Cx{}] = 9;

  data2::SoA<Pos const&, Cx const&> soa;

  data2::SoA<Pos&, Cx&> soa2{test};

  soa2[Pos{}] = p[Pos{}];

  //   Eigen::internal::set_is_malloc_allowed(true);
  soa = soa2;
  //   Eigen::internal::set_is_malloc_allowed(false);

  pretty_print("soa p", soa[Pos{}]);
  pretty_print("soa x", soa[Cx{}]);

  //   return 0;

  return 0;

  constexpr size_t N = 6;

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