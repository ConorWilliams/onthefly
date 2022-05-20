

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <iostream>
#include <random>

#include "libatom/io/xyz.hpp"
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
    atoms(AtomicNum{}, i) = 26;
    atoms(Frozen{}, i) = false;
  }

  return atoms;
}

#include "libatom/neighbour/lcr_sort.hpp"

auto main(int, char**) -> int {
  //

  SimCell atoms({{34, 34, 34}, {true, true, true}});

  random_simcell(atoms, 10'000);

  static_cast<AtomArray<Position, Frozen, AtomicNum>>(atoms)
      = lcr_sort({atoms.box, 6, false}, atoms);

  fmt::print("num atoms is {}\n", atoms.size());

  floating rcut = 6;

  auto out = fmt::output_file("dump.xyz");

  dump_xyz(out, atoms, "My comment");

  {
    NeighbourList neigh(atoms.box, rcut);

    neigh.rebuild(atoms);  // Warm up + alloc

    timeit("Fast", [&] { neigh.rebuild_parallel(atoms); });

    std::cout << "working\n";

    return 0;
  }
}