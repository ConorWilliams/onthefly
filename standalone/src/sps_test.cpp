#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

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

  cell(AtomicNum{}, cell.size() - 1) = 1;  // Set last is hydrogen

  return cell;
}

auto main(int, char **) -> int {
  //

  std::ifstream file("/home/cdt1902/phd/onthefly/data/xyz/V1.xyz", std::ios::in);

  SimCell cell = init_cell(431, 17.1598406);

  env::EnvCell envs({5.2}, cell);

  io::stream_xyz(file, cell);

  io::stream_xyz(file, cell);

  //   {  //
  //     env::EnvCell classifyer({5.2}, atoms);

  //     classifyer.rebuild(atoms, omp_get_max_threads());
  //   }

  //   {
  //     // std::cout << atoms(Axis{}) << std::endl;

  //     potentials::EAM
  //     pot{std::make_shared<potentials::DataEAM>(std::ifstream{"../data/wen.eam.fs"})};

  //     neighbour::List nl(atoms, pot.rcut());

  //     nl.rebuild(atoms, omp_get_max_threads());

  //     floating E0 = pot.energy(atoms, nl, omp_get_max_threads());

  //     fmt::print("E0={}\n", E0);

  //     saddle::perturb({2.8, 2.8, 2.8}, atoms, 4, 0.6);

  //     // saddle::perturb({1.54758, 1.48733, 4.30367}, atoms, 4, 0.6);

  //     potentials::Dimer::Options A;

  //     A.debug = true;
  //     // A.iter_max_rot = 100;
  //     // A.relax_in_convex = false;
  //     // A.theta_tol = 5 * 2 * 3.141 / 360;

  //     potentials::Dimer dim(A, std::make_unique<potentials::EAM>(pot));

  //     minimise::LBFGS::Options opt;

  //     opt.debug = true;
  //     // opt.convex_max = 7;

  //     minimise::LBFGS lbfgs(opt);

  //     if (lbfgs.minimise(atoms, dim, omp_get_max_threads())) {
  //       io::dump_xyz(fmt::output_file("saddle.xyz"), atoms, "");

  //       nl.rebuild(atoms, omp_get_max_threads());

  //       floating Ef = pot.energy(atoms, nl, omp_get_max_threads());

  //       fmt::print("Ef-E0={}\n", Ef - E0);
  //     }
  //   }

  return 0;
}