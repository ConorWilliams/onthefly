

#include <fmt/core.h>

#include <iostream>
#include <random>

#include "libatom/neighbour/neighbour_list.hpp"
#include "libatom/system/atom_array.hpp"
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

auto main(int, char**) -> int {
  //

  SimCell atoms({{17, 17, 17}, {true, true, true}});

  random_simcell(atoms, 1000);

  fmt::print("num atoms is {}\n", atoms.size());

  floating rcut = 6;

  {
    NeighbourCell neigh;

    neigh.rebuild_neighbour_lists(atoms, rcut);  // Warm up + alloc

    timeit("Fast", [&] { neigh.rebuild_neighbour_lists(atoms, rcut); });

    std::cout << "working\n";

    return 0;
  }
}