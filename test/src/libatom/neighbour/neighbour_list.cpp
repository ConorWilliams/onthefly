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
  NeighbourCell neigh;

  neigh.rebuild_neighbour_lists(atoms, rcut);

  std::vector<std::vector<Neigh>> nl;

  slow_neigh_list(nl, atoms, rcut);

  for (size_t i = 0; i < atoms.size(); i++) {
    //
    std::vector<Neigh> nl2;

    neigh.for_neighbours(i, rcut, [&](std::size_t n, floating, Vec3<floating> const& dr) {
      nl2.push_back({neigh.image_to_real(n), dr});
    });

    std::sort(nl2.begin(), nl2.end(), [](auto const& a, auto const& b) { return a.i < b.i; });

    for (size_t j = 0; j < nl2.size(); j++) {
      // Test same number of neighbours
      REQUIRE(nl2.size() == nl[i].size());
      // Test all neighbours have the same index
      REQUIRE(nl2[j].i == nl[i][j].i);
      // Test all neighbours have the same minimum image positions
      REQUIRE(norm(nl2[j].dr - nl[i][j].dr) < 0.001);
    }
  }
}

TEST_CASE("Neighbour list fuzz testing") {
  for (size_t i = 0; i < 1000; i++) {
    //
    SimCell atoms({10 * dis(gen), 10 * dis(gen), 10 * dis(gen)},
                  {dis(gen) < 0.5, dis(gen) < 0.5, dis(gen) < 0.5});
  }
}