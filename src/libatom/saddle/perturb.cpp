#include <random>

#include "libatom/asserts.hpp"
#include "libatom/random/xoshiro.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::saddle {

void perturb(Vec3<floating> const& centre, SimCell& cell, floating rcut, floating stddev) {
  // PRNG
  static thread_local std::random_device rd;
  static thread_local random::Xoshiro rng({rd(), rd(), rd(), rd()});
  // static thread_local random::Xoshiro rng({1, 1, 1, 1});

  std::normal_distribution<floating> normal(0, 1);
  std::normal_distribution<floating> gauss(0, stddev);

  cell(Axis{}) = 0; // Zero the rotatio axis.

  for (std::size_t i = 0; i < cell.size(); i++) {
    //
    if (!cell(Frozen{}, i)) {
      floating dr2 = norm_sq(cell.min_image(cell(Position{}, i), centre));

      if (dr2 < rcut * rcut) {
        //
        floating env = 1 - std::sqrt(dr2) / rcut;

        cell(Position{}, i) += env * Vec3<floating>{gauss(rng), gauss(rng), gauss(rng)};
        cell(Axis{}, i) += Vec3<floating>{normal(rng), normal(rng), normal(rng)};
      }
    }
  }

  cell(Axis{}) /= norm(cell(Axis{}));
}

} // namespace otf::saddle