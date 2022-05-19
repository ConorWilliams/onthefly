#include "libatom/neighbour/neighbour_list.hpp"

#include <algorithm>
#include <cstddef>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "fmt/core.h"
#include "libatom/system/member.hpp"
#include "libatom/system/sim_cell.hpp"
#include "libatom/utils.hpp"

using namespace otf;

std::mt19937 gen(33);
std::uniform_real_distribution<floating> dis(0, 1);

SimCell random_simcell(SimCell& atoms, std::size_t n) {
  //
  atoms.resize(n);

  for (size_t i = 0; i < n; i++) {
    atoms(Position{}, i) = Vec3<floating>{dis(gen), dis(gen), dis(gen)} * atoms.box.extents();
    atoms(AtomicNum{}, i) = 1;
    atoms(Frozen{}, i) = false;
  }

  return atoms;
}

struct Neigh {
  std::size_t i;
  Vec3<floating> dr;
};

void slow_neigh_list(std::vector<std::vector<Neigh>>& nl, SimCell const& atoms, floating rcut) {
  nl.resize(atoms.size());

  for (size_t i = 0; i < atoms.size(); i++) {
    nl[i].clear();

    for (size_t j = 0; j < atoms.size(); j++) {
      if (i != j) {
        Vec3<floating> dr = atoms.box.min_image(atoms(Position{}, i), atoms(Position{}, j));

        if (norm_sq(dr) < rcut * rcut) {
          nl[i].push_back({j, dr});
        }
      }
    }
    std::sort(nl[i].begin(), nl[i].end(), [](auto const& a, auto const& b) { return a.i < b.i; });
  }
}

void test(SimCell const& atoms, floating rcut) {
  //
  NeighbourList neigh;

  neigh.rebuild(atoms, rcut);

  std::vector<std::vector<Neigh>> nl;

  slow_neigh_list(nl, atoms, rcut);

  for (size_t i = 0; i < atoms.size(); i++) {
    //
    std::vector<Neigh> nl2;

    neigh.for_neighbours(i, [&](std::size_t n, floating, Vec3<floating> const& dr) {
      nl2.push_back({neigh.image_to_real(n), dr});
    });

    std::sort(nl2.begin(), nl2.end(), [](auto const& a, auto const& b) { return a.i < b.i; });

    // Test same number of neighbours
    REQUIRE(nl2.size() == nl[i].size());

    for (size_t j = 0; j < nl2.size(); j++) {
      // Test all neighbours have the same index
      REQUIRE(nl2[j].i == nl[i][j].i);
      // Test all neighbours have the same minimum image positions
      REQUIRE(norm(nl2[j].dr - nl[i][j].dr) < 0.001);
    }
  }
}

TEST_CASE("Neighbour list speed testing") {
  SimCell atoms({{1, 2, 1}, {true, true, true}});

  random_simcell(atoms, 10'000 * atoms.box.extents().prod());

  fmt::print("num atoms is {}\n", atoms.size());

  floating rcut = 0.1;

  {
    NeighbourList neigh;

    neigh.rebuild(atoms, rcut);  // Warm up + alloc

    timeit("Fast", [&] { neigh.rebuild(atoms, rcut); });

    int x = 0;
    int y = 0;

    timeit("Counting", [&] {
      for (size_t i = 0; i < atoms.size(); i++) {
        //
        x++;
        neigh.for_neighbours(i, [&](std::size_t, floating, Vec3<floating> const&) { y++; });
      }
    });

    fmt::print("num neigh = {}\n", (double)y / x);
  }

  {
    std::vector<std::vector<Neigh>> nl;

    slow_neigh_list(nl, atoms, rcut);  // Warm up + alloc

    timeit("Slow", [&] { slow_neigh_list(nl, atoms, rcut); });
  }
}

TEST_CASE("Neighbour list fuzz testing") {
  for (size_t i = 0; i < 100; i++) {
    //
    SimCell atoms({
        {1 + 9 * dis(gen), 1 + 9 * dis(gen), 1 + 9 * dis(gen)},
        {dis(gen) < 0.5, dis(gen) < 0.5, dis(gen) < 0.5},
    });

    random_simcell(atoms, 100 + 900 * dis(gen));

    for (size_t j = 0; j < 10; j++) {
      test(atoms, 0.25 + 0.2499 * dis(gen));
    }
  }
}