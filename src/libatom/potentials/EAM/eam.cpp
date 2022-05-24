

#include "libatom/potentials/EAM/eam.hpp"

#include <cmath>
#include <cstddef>
#include <memory>

#include "libatom/asserts.hpp"
#include "libatom/neighbour/neighbour_list.hpp"
#include "libatom/potentials/EAM/data.hpp"
#include "libatom/system/atom_array.hpp"
#include "libatom/system/member.hpp"
#include "libatom/system/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf {

  floating EAM::energy(SimCell const &x, NeighbourList const &nl, std::size_t num_threads) const {
    floating v_sum = 0;
    floating f_sum = 0;

#pragma omp parallel for reduction(+ : v_sum, f_sum) num_threads(num_threads) schedule(static)
    for (size_t a = 0; a < x.size(); a++) {
      // Skip contibution from frozen atoms
      if (x(Frozen{}, a)) {
        continue;
      }

      floating rho = 0;

      nl.for_neighbours(a, rcut(), [&](auto bp, floating r_sq, Vec3<floating> const &) {
        //
        floating r = std::sqrt(r_sq);
        std::size_t b = nl.image_to_real(bp);

        v_sum += m_data->v(x(AtomicNum{}, a), x(AtomicNum{}, b)).f(r);
        rho += m_data->phi(x(AtomicNum{}, b), x(AtomicNum{}, a)).f(r);
      });

      f_sum += m_data->f(x(AtomicNum{}, a)).f(rho);
    }

    return (0.5 * v_sum) + f_sum;
  }

  void EAM::gradient(SimCell &x, NeighbourList const &nl, std::size_t num_threads) {
    // Usually a noop
    m_aux.destructive_resize(x.size());

// First sum computes density  at each atom, runs over active+boundary atoms
#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (size_t b = 0; b < x.size(); b++) {
      floating rho = 0;

      // Computes rho at atom
      nl.for_neighbours(b, [&](std::size_t a, floating r_sq, Vec3<floating> const &) {
        rho += m_data->phi(x(AtomicNum{}, nl.image_to_real(a)), x(AtomicNum{}, b))
                   .f(std::sqrt(r_sq));
      });

      // Compute F'(rho) at atom
      m_aux(Fprime{}, b) = m_data->f(x(AtomicNum{}, b)).fp(rho);
    }

// Second sum computes gradient, only runs over active atoms
#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (std::size_t g = 0; g < x.size(); ++g) {
      //
      Vec3<floating> grad = Vec3<floating>::Zero();

      if (!x(Frozen{}, g)) {
        nl.for_neighbours(g, [&](std::size_t ap, floating r_sq, Vec3<floating> const &dr) {
          //
          floating r = std::sqrt(r_sq);
          std::size_t a = nl.image_to_real(ap);

          floating mag
              = m_data->v(x(AtomicNum{}, a), x(AtomicNum{}, g)).fp(r)
                + m_aux(Fprime{}, g) * m_data->phi(x(AtomicNum{}, a), x(AtomicNum{}, g)).fp(r)
                + m_aux(Fprime{}, a) * m_data->phi(x(AtomicNum{}, g), x(AtomicNum{}, a)).fp(r);

          grad -= (mag / r) * dr;
        });
      }

      // Write grad to atom
      x(Gradient{}, g) = grad;
    }
  }

  void EAM::hessian(SimCell const &, NeighbourList const &, std::size_t) const {}

}  // namespace otf