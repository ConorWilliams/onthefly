#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "libatom/atom.hpp"
#include "libatom/env/geometry.hpp"
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

  Vec shape{6, 6, 6};  // In unit cells

  for (std::size_t k = 0; k < shape[0]; k++) {
    for (std::size_t j = 0; j < shape[1]; j++) {
      for (std::size_t i = 0; i < shape[2]; i++) {
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

  for (std::size_t i = 0; i < lat.size(); i++) {
    atoms(Position{}, i) = lat[i].off;
    atoms(AtomicNum{}, i) = lat[i].num;
    atoms(Frozen{}, i) = false;
  }

  //   atoms(Frozen{}, 341) = true;
  //   atoms(Frozen{}, 0) = true;

  auto f = fmt::output_file("lbfgs_debug.xyz");

  io::dump_xyz(f, atoms, fmt::format("Temp={}", T));

  fmt::print("Frozen = {}\n", atoms.count_frozen());

  AtomVector<Position, AtomicNum> x;  // Initialise a vector of 0 atoms.

  x.emplace_back({1, 2, 3}, 11);  // Ad an atom

  Vec3<floating> p = x[0](Position{});

  fmt::print("{},{}\n", x[0](AtomicNum{}), p[0]);

  // exit(1);

  {
    // saddle::perturb({7, 7, 7}, atoms, 4, 0.6);
    saddle::perturb({2.8, 2.8, 2.8}, atoms, 4, 0.6);

    // std::cout << atoms(Axis{}) << std::endl;

    potentials::EAM pot{std::make_shared<potentials::DataEAM>(std::ifstream{"../data/wen.eam.fs"})};

    std::unique_ptr<potentials::Base> tmp = std::make_unique<potentials::EAM>(std::move(pot));

    potentials::Dimer::Options A;

    A.debug = true;
    // A.iter_max_rot = 100;
    // A.relax_in_convex = false;
    // A.theta_tol = 5 * 2 * 3.141 / 360;

    potentials::Dimer dim(A, std::move(tmp));

    minimise::LBFGS::Options opt;

    opt.debug = true;
    // opt.convex_max = 7;

    minimise::LBFGS lbfgs(opt);

    if (lbfgs.minimise(atoms, dim, omp_get_max_threads())) {
      io::dump_xyz(fmt::output_file("saddle.xyz"), atoms, "");
    }
  }

  return 0;
}