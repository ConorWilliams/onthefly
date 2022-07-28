

#include "libatom/potentials/EAM/eam.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>

#include "fmt/core.h"
#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/neighbour/list.hpp"
#include "libatom/potentials/EAM/data.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::potentials {

floating EAM::energy(SimCell const& x, neighbour::List& nl, std::size_t num_threads) {
  floating v_sum = 0;
  floating f_sum = 0;

#pragma omp parallel for reduction(+ : v_sum, f_sum) num_threads(num_threads) schedule(static)
  for (std::size_t a = 0; a < x.size(); a++) {
    // Skip contibution from frozen atoms
    if (x(Frozen{}, a)) {
      continue;
    }

    floating rho = 0;

    nl.for_neighbours(a, rcut(), [&](auto bp, floating r, Vec3<floating> const&) {
      //
      std::size_t b = nl.image_to_real(bp);

      v_sum += m_data->v(x(AtomicNum{}, a), x(AtomicNum{}, b)).f(r);
      rho += m_data->phi(x(AtomicNum{}, b), x(AtomicNum{}, a)).f(r);
    });

    f_sum += m_data->f(x(AtomicNum{}, a)).f(rho);
  }

  return (0.5 * v_sum) + f_sum;
}

std::optional<floating> EAM::gradient(SimCell& x, neighbour::List& nl, std::size_t num_threads) {
  // Usually a noop
  m_aux.destructive_resize(x.size());

// First sum computes density  at each atom, runs over active+boundary atoms
#pragma omp parallel for num_threads(num_threads) schedule(static)
  for (std::size_t b = 0; b < x.size(); b++) {
    floating rho = 0;

    // Computes rho at atom
    nl.for_neighbours(b, rcut(), [&](std::size_t a, floating r, Vec3<floating> const&) {
      rho += m_data->phi(x(AtomicNum{}, nl.image_to_real(a)), x(AtomicNum{}, b)).f(r);
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
      nl.for_neighbours(g, rcut(), [&](std::size_t ap, floating r, Vec3<floating> const& dr) {
        //
        std::size_t a = nl.image_to_real(ap);

        floating mag = m_data->v(x(AtomicNum{}, a), x(AtomicNum{}, g)).fp(r) +
                       m_aux(Fprime{}, g) * m_data->phi(x(AtomicNum{}, a), x(AtomicNum{}, g)).fp(r) +
                       m_aux(Fprime{}, a) * m_data->phi(x(AtomicNum{}, g), x(AtomicNum{}, a)).fp(r);

        grad -= (mag / r) * dr;
      });
    }

    // Write grad to atom
    x(Gradient{}, g) = grad;
  }

  return std::nullopt;
}

void EAM::hessian(SimCell& x, neighbour::List& nl, std::size_t) {
  // Usually a noop, make space in aux
  m_aux.destructive_resize(x.size());

  // Compute hessian indexes
  std::size_t count = 0;
  for (std::size_t i = 0; i < x.size(); i++) {
    if (!x(Frozen{}, i)) {
      m_aux(Hidx{}, i) = count++;
    }
  }

  // First sum computes  rho & mu  at each atom, runs over active + frozen atoms
  for (std::size_t b = 0; b < x.size(); ++b) {
    floating rho = 0;
    Vec3<floating> mu = Vec3<floating>::Zero();
    // Compute rho and mu at atom via local sum
    nl.for_neighbours(b, rcut(), [&](std::size_t ap, floating r, Vec3<floating> const& dr) {
      //
      std::size_t a = nl.image_to_real(ap);

      rho += m_data->phi(x(AtomicNum{}, a), x(AtomicNum{}, b)).f(r);
      mu -= m_data->phi(x(AtomicNum{}, a), x(AtomicNum{}, b)).fp(r) * dr / r; // A.13
    });

    // Write
    m_aux(Rho{}, b) = rho;
    m_aux(Mu{}, b) = mu;
  }

  x.zero_hess();

  // Second sums computes hessian, running over all atoms
  for (std::size_t z = 0; z < x.size(); ++z) {
    // During this section we only write to the z^th column triplet

    if (!x(Frozen{}, z)) {
      auto v = m_aux(Mu{}, z).matrix();
      // A.15 pre sum term
      x.hess(m_aux(Hidx{}, z), m_aux(Hidx{}, z)).matrix() = m_data->f(x(AtomicNum{}, z)).fpp(m_aux(Rho{}, z)) * v * v.transpose();

      // Note: dr = r^{za}
      nl.for_neighbours(z, rcut(), [&](std::size_t ap, floating r, Vec3<floating> const& dr) {
        // First compute sum over neigh for block-diagonal hessian-elements, neighbours can be
        // frozen or not

        std::size_t a = nl.image_to_real(ap);

        // A.14
        floating A = (m_data->v(x(AtomicNum{}, a), x(AtomicNum{}, z)).fp(r) +

                      m_data->f(x(AtomicNum{}, z)).fp(m_aux(Rho{}, z)) * m_data->phi(x(AtomicNum{}, a), x(AtomicNum{}, z)).fp(r) +

                      m_data->f(x(AtomicNum{}, a)).fp(m_aux(Rho{}, a)) * m_data->phi(x(AtomicNum{}, z), x(AtomicNum{}, a)).fp(r)) /
                     r;

        // A.14
        floating B = m_data->v(x(AtomicNum{}, a), x(AtomicNum{}, z)).fpp(r) +

                     m_data->f(x(AtomicNum{}, z)).fp(m_aux(Rho{}, z)) * m_data->phi(x(AtomicNum{}, a), x(AtomicNum{}, z)).fpp(r) +

                     m_data->f(x(AtomicNum{}, a)).fp(m_aux(Rho{}, a)) * m_data->phi(x(AtomicNum{}, z), x(AtomicNum{}, a)).fpp(r);

        floating ABFpp =
            (A - B -
             m_data->f(x(AtomicNum{}, a)).fpp(m_aux(Rho{}, a)) * ipow<2>(m_data->phi(x(AtomicNum{}, z), x(AtomicNum{}, a)).fp(r))) /
            (r * r);

        x.hess(m_aux(Hidx{}, z), m_aux(Hidx{}, z)).matrix() +=
            A * Eigen::Matrix<floating, spatial_dims, spatial_dims>::Identity() - ABFpp * dr.matrix() * dr.matrix().transpose();

        // Now we will compute the off diagonal element in this column block of H, the a's must
        // now not be frozen, we do not commpute overlap here
        if (!x(Frozen{}, a)) {
          floating BArr = (B - A) / (r * r);

          floating ur =
              m_data->f(x(AtomicNum{}, z)).fpp(m_aux(Rho{}, z)) * m_data->phi(x(AtomicNum{}, a), x(AtomicNum{}, z)).fp(r) / r;

          floating ru =
              m_data->f(x(AtomicNum{}, a)).fpp(m_aux(Rho{}, a)) * m_data->phi(x(AtomicNum{}, z), x(AtomicNum{}, a)).fp(r) / r;

          // TODO : .noaliase()

          x.hess(m_aux(Hidx{}, z), m_aux(Hidx{}, a)).matrix() =
              -BArr * dr.matrix() * dr.matrix().transpose() - A * Eigen::Matrix<floating, spatial_dims, spatial_dims>::Identity() +
              ur * m_aux(Mu{}, z).matrix() * dr.matrix().transpose() - ru * m_aux(Mu{}, a).matrix() * dr.matrix().transpose();
        }
      });
    }
  }

  // Finally we must compute the overlap terms
  for (std::size_t a = 0; a < x.size(); ++a) {
    //
    floating ddFg = m_data->f(x(AtomicNum{}, a)).fpp(m_aux(Rho{}, a));

    nl.for_neighbours(a, rcut(), [&](std::size_t ep, floating r_e, auto const& dr_ae) {
      nl.for_neighbours(a, rcut(), [&](std::size_t gp, floating r_g, auto const& dr_ag) {
        std::size_t e = nl.image_to_real(ep);
        std::size_t g = nl.image_to_real(gp);
        if (e != g && !x(Frozen{}, e) && !x(Frozen{}, g)) {
          // Now iterating over all pair of unfrozen neighbours of a
          floating mag = ddFg / (r_e * r_g) * m_data->phi(x(AtomicNum{}, g), x(AtomicNum{}, a)).fp(r_g) *
                         m_data->phi(x(AtomicNum{}, e), x(AtomicNum{}, a)).fp(r_e);

          x.hess(m_aux(Hidx{}, e), m_aux(Hidx{}, g)).matrix() += mag * dr_ae.matrix() * dr_ag.matrix().transpose();
        }
      });
    });
  }
}
} // namespace otf::potentials