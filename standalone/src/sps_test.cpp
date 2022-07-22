#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/env/geometry.hpp"
#include "libatom/env/topology.hpp"
#include "libatom/io/xyz.hpp"
#include "libatom/minimise/LBFGS/core.hpp"
#include "libatom/minimise/LBFGS/lbfgs.hpp"
#include "libatom/neighbour/list.hpp"
#include "libatom/potentials/EAM/eam.hpp"
#include "libatom/potentials/ROT/dimer.hpp"
#include "libatom/random/xoshiro.hpp"
#include "libatom/saddle/perturb.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

using namespace otf;

SimCell init_cell(int n, floating l) {
  //
  SimCell cell({{l, l, l}, {true, true, true}});

  cell.destructive_resize(n);

  for (std::size_t i = 0; i < cell.size(); i++) {
    cell(Frozen{}, i) = false;
    cell(AtomicNum{}, i) = 26;
  }

  std::ifstream file("/home/cdt1902/phd/onthefly/data/xyz/V1.xyz", std::ios::in);

  io::stream_xyz(file, cell);

  return cell;
}

auto main(int, char **) -> int {
  //

  SimCell atoms = init_cell(431, 17.1598406);

  potentials::EAM pot{std::make_shared<potentials::DataEAM>(std::ifstream{"../data/FeH-A.fs"})};

  minimise::LBFGS lbfgs = [&]() {
    //
    minimise::LBFGS::Options opt;

    opt.debug = true;

    return minimise::LBFGS{opt};
  }();

  neighbour::List nl(atoms, pot.rcut());

  floating E0 = [&] {
    //
    fmt::print("Initial minimisation...\n");

    bool x = time_call("min", &minimise::LBFGS::minimise, lbfgs, atoms, pot, omp_get_max_threads());

    ASSERT(x, "initial min failed");

    nl.rebuild(atoms, omp_get_max_threads());

    return pot.energy(atoms, nl, omp_get_max_threads());
  }();

  fmt::print("E0={}\n", E0);

  floating Ef = [&] {
    //
    potentials::Dimer::Options opt;

    opt.debug = true;

    potentials::Dimer dim{opt, std::make_unique<potentials::EAM>(pot)};

    while (true) {
      auto copy = atoms;

      //   saddle::perturb({2.8, 2.8, 2.8}, copy, 4, 0.6);

      saddle::perturb({3.85852, 3.85841, 6.75341}, copy, 4, 0.6);

      if (time_call("SPS", &minimise::LBFGS::minimise, lbfgs, copy, dim, omp_get_max_threads())) {
        //
        nl.rebuild(copy, omp_get_max_threads());

        return pot.energy(copy, nl, omp_get_max_threads());

      } else {
        fmt::print("SPS-fail\n");
      }
    }
  }();

  fmt::print("Ef-E0={}\n", Ef - E0);

  return 0;
}