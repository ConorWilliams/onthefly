

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <random>

#include "libatom/atom_array.hpp"
#include "libatom/io/xyz.hpp"
#include "libatom/neighbour/list.hpp"
#include "libatom/potentials/EAM/eam.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

using namespace otf;

std::mt19937 gen(33);
std::uniform_real_distribution<floating> dis(0, 1);

SimCell random_simcell(SimCell& atoms, std::size_t n) {
  //
  atoms.destructive_resize(n);

  for (std::size_t i = 0; i < n; i++) {
    atoms(Position{}, i) = Vec3<floating>{dis(gen), dis(gen), dis(gen)} * atoms.extents();
    atoms(AtomicNum{}, i) = 26;
    atoms(Frozen{}, i) = false;
  }

  return atoms;
}

auto main(int, char**) -> int {
  //

  SimCell atoms({{17, 17, 17}, {true, true, true}});

  random_simcell(atoms, 1'000);

  fmt::print("num atoms is {}\n", atoms.size());

  neighbour::List neigh(atoms, 6);

  neigh.rebuild(atoms, omp_get_max_threads());

  timeit("Fast", [&] { neigh.rebuild(atoms, omp_get_max_threads()); });

  return 0;
}