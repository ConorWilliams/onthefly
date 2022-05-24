#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "libatom/io/xyz.hpp"
#include "libatom/minimise/LBFGS/core.hpp"
#include "libatom/minimise/LBFGS/lbfgs.hpp"
#include "libatom/neighbour/neighbour_list.hpp"
#include "libatom/potentials/EAM/eam.hpp"
#include "libatom/system/member.hpp"
#include "libatom/system/sim_cell.hpp"
#include "libatom/utils.hpp"

using namespace otf;

auto main(int, char **) -> int {
  //
  using Vec = Vec3<floating>;

  struct MotifPt {
    std::size_t num;
    Vec off;
  };

  // Fractional motif
  std::array BCC = {
      MotifPt{26, Vec{0.0, 0.0, 0.0}},
      MotifPt{26, Vec{0.5, 0.5, 0.5}},
  };

  floating T = 300.0;

  floating V
      = 11.64012 + T * (9.37798e-5 + T * (3.643134e-7 + T * (1.851593e-10 + T * 5.669148e-14)));

  floating a = std::pow(2 * V, 1.0 / 3);  // lat param

  fmt::print("lattice parameter = {}\n", a);

  std::vector<MotifPt> lat;

  Vec shape{10, 10, 10};  // In unit cells

  for (size_t k = 0; k < shape[0]; k++) {
    for (size_t j = 0; j < shape[1]; j++) {
      for (size_t i = 0; i < shape[2]; i++) {
        for (auto &&mot : BCC) {
          //
          Vec lp = Vec{(double)i, (double)j, (double)k} + mot.off;

          lat.push_back({mot.num, lp * a});
        }
      }
    }
  }

  //   std::random_shuffle(lat.begin(), lat.end());

  lat.erase(lat.begin() + 1);

  SimCell atoms{{a * shape, {true, true, true}}};

  atoms.destructive_resize(lat.size());

  for (size_t i = 0; i < lat.size(); i++) {
    atoms(Position{}, i) = lat[i].off;
    atoms(AtomicNum{}, i) = lat[i].num;
    atoms(Frozen{}, i) = false;
  }

  atoms(Frozen{}, 0) = true;

  auto f = fmt::output_file("dump.xyz");

  dump_xyz(f, atoms, fmt::format("Temp={}", T));

  {
    EAM pot{std::make_shared<DataEAM>(std::ifstream{"../data/wen.eam.fs"})};

    NeighbourList neigh(atoms.box, pot.rcut() + 1);

    neigh.rebuild(atoms, omp_get_max_threads());

    double energy = 0;

    fmt::print("num threads = {}\n", omp_get_max_threads());

    timeit("Energy call", [&] { energy = pot.energy(atoms, neigh, omp_get_max_threads()); });

    fmt::print("Energy = {}\n", energy);

    timeit("Grad call", [&] { pot.gradient(atoms, neigh, omp_get_max_threads()); });

    LBFGS::Options opt;

    opt.debug = false;
    opt.skin_frac = 1.1;

    LBFGS lbfgs(opt);

    timeit("Minimise", [&] {
      auto copy = atoms;
      lbfgs.minimise(copy, pot, omp_get_max_threads());
    });
  }

  return 0;
}