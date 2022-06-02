#pragma once

#include <cmath>
#include <cstddef>
#include <memory>

#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/minimise/LBFGS/core.hpp"
#include "libatom/neighbour/list.hpp"
#include "libatom/potentials/base.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::potentials {

  /**
   * @brief Dimer is a potential adaptor.
   *
   * It wraps a potential and inverts the commponent of the gradient parallel to the minimum mode.
   */
  class Dimer : public Base {
  public:
    /**
     * @brief Used to configure the dimers internal minimisation pass.
     */
    struct Options {
      /** @brief L-BFGS core rotation-history size */
      std::size_t n = 6;
      /** @brief Maximum number of iterations during rotation minimization */
      std::size_t iter_max_rot = 20;
      /** @brief Half dimer length */
      double delta_r = 0.001;
      /** @brief (Rad) rotation convergence criterion */
      double theta_tol = 2 * M_PI / 360.;
      /** @brief If true when convex we return only the component parallel to the min mode. */
      bool relax_in_convex = true;
      /** @brief Print out debug info */
      bool debug = false;
    };

    Dimer(Options opt, std::unique_ptr<Base> to_wrap)
        : m_opt(opt), m_core(opt.n), m_wrapped(std::move(to_wrap)) {}

    /**
     * @brief Copies this potentials cut-off radius
     */
    std::unique_ptr<Base> clone() const override {
      return std::make_unique<Dimer>(m_opt, m_wrapped->clone());
    }

    /**
     * @brief Get this potentials cut-off radius that the neighbour lists should be configured for.
     *
     * Here we retun double the wrapped cut-off plus two times delta_r, this means we can displace
     * every atom by delta_r and the not need to rebuild the neighbour lists.
     */
    floating rcut() const override { return m_wrapped->rcut() + m_opt.delta_r * 2; }

    /**
     * @brief Dimer does not support energy.
     */
    [[noreturn]] floating energy(SimCell const &, neighbour::List &, std::size_t) override {
      throw unsupported{};
    }

    /**
     * @brief Compute gradient, assumes the neighbour list are ready.
     *
     * Force on frozen atoms will be zero.
     */
    std::optional<floating> gradient(SimCell &, neighbour::List &, std::size_t num_thread) override;

    /**
     * @brief Dimer does not support hessians.
     */
    [[noreturn]] void hessian(SimCell &, neighbour::List &, std::size_t) override {
      throw unsupported{};
    }

  private:
    Options m_opt;
    minimise::CoreLBFGS m_core;
    std::unique_ptr<Base> m_wrapped;

    Position::matrix_type m_delta;       // Store displacement for updates
    Position::matrix_type m_delta_prev;  // Store previous displacement for updates
    Position::matrix_type m_axisp;       // Temporary axis

    Gradient::matrix_type m_g0;       // Central grad
    Gradient::matrix_type m_g1;       // Temp end grad
    Gradient::matrix_type m_delta_g;  // Gradient difference
  };

}  // namespace otf::potentials