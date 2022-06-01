

#include "libatom/potentials/ROT/dimer.hpp"

#include <fmt/core.h>

#include <cmath>
#include <cstddef>
#include <memory>
#include <type_traits>

#include "Eigen/src/Core/util/Constants.h"
#include "libatom/asserts.hpp"
#include "libatom/atom_array.hpp"
#include "libatom/minimise/LBFGS/core.hpp"
#include "libatom/neighbour/list.hpp"
#include "libatom/potentials/base.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::potentials {

  std::optional<floating> Dimer::gradient(SimCell &cell, neighbour::List &nl,
                                          std::size_t num_threads) {
    //
    m_core.clear();

    m_wrapped->gradient(cell, nl, num_threads);  // Gradient at centre (g1)
    m_g0 = cell(Gradient{});

    m_delta = -m_opt.delta_r * cell(Axis{});
    nl.update_positions(m_delta);
    using std::swap;  // ADL
    swap(m_delta, m_delta_prev);

    m_wrapped->gradient(cell, nl, num_threads);  // Gradient at end (g1)
    m_g1 = cell(Gradient{});

    floating curv = [&] {
      for (std::size_t i = 0;; i++) {
        m_delta_g = m_g1 - m_g0;
        m_delta_g -= gdot(m_delta_g, cell(Axis{})) * cell(Axis{});  // Torque

        // Use lbfgs to find rotation plane
        auto &theta = m_core.newton_step(cell(Axis{}), m_delta_g);

        theta -= gdot(theta, cell(Axis{})) * cell(Axis{});  // Ortho
        theta *= 1 / norm(theta);                           //      normalization

        floating b_1 = gdot(m_g1 - m_g0, theta) / m_opt.delta_r;
        floating c_x0 = gdot(m_g1 - m_g0, cell(Axis{})) / m_opt.delta_r;
        floating theta_1 = -0.5 * std::atan(b_1 / std::abs(c_x0));  // Trial rotation angle

        if (std::abs(theta_1) < m_opt.theta_tol || i == m_opt.iter_max_rot) {
          return c_x0;
        } else {
          // Trial rotation
          m_axisp = cell(Axis{}) * std::cos(theta_1) + theta * std::sin(theta_1);

          m_delta = -m_opt.delta_r * m_axisp;  // Temporarily store next m_delta_prev into m_delta
          swap(m_delta, m_delta_prev);         // Now put it in the correct place
          m_delta = m_delta_prev - m_delta;
          nl.update_positions(m_delta);

          m_wrapped->gradient(cell, nl, num_threads);  // Gradient at primed end (g1p)

          double c_x1 = gdot(cell(Gradient{}) - m_g0, m_axisp) / m_opt.delta_r;
          double a_1 = (c_x0 - c_x1 + b_1 * sin(2 * theta_1)) / (1 - std::cos(2 * theta_1));
          double theta_min = 0.5 * std::atan(b_1 / a_1);  // Optimal rotation

          // Flip if extrema is maxima
          if (a_1 * std::cos(2 * theta_min) - a_1 + b_1 * std::sin(2 * theta_min) > 0) {
            theta_min += M_PI / 2;
          }

          cell(Axis{}) = cell(Axis{}) * std::cos(theta_min) + theta * std::sin(theta_min);

          // Interpolate force at new rotation
          m_g1 = (std::sin(theta_1 - theta_min) / std::sin(theta_1)) * m_g1
                 + (std::sin(theta_min) / std::sin(theta_1)) * cell(Gradient{})
                 + (1 - std::cos(theta_min) - std::sin(theta_min) * std::tan(0.5 * theta_1)) * m_g0;

          if (m_opt.debug) {
            constexpr auto str = "    Dimer: i={:<4} theta={:f} curv={:f}\n";
            fmt::print(str, i, std::abs(theta_min),
                       gdot(m_g1 - m_g0, cell(Axis{})) / m_opt.delta_r);
          }

          if (std::abs(theta_min) < m_opt.theta_tol) {
            return gdot(m_g1 - m_g0, cell(Axis{})) / m_opt.delta_r;
          }
        }
      }
    }();

    m_delta = -m_delta_prev;
    nl.update_positions(m_delta);

    if (!m_opt.relax_in_convex && curv > 0) {
      cell(Gradient{}) = -gdot(m_g0, cell(Axis{})) * cell(Axis{});
    } else {
      cell(Gradient{}) = m_g0 - 2 * gdot(m_g0, cell(Axis{})) * cell(Axis{});
    }

    return curv;
  }

}  // namespace otf::potentials