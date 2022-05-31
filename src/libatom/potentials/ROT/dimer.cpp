

#include "libatom/potentials/ROT/dimer.hpp"

#include <cstddef>
#include <memory>
#include <type_traits>
#include <cmath>

#include "libatom/asserts.hpp"
#include "libatom/atom_array.hpp"
#include "libatom/minimise/LBFGS/core.hpp"
#include "libatom/neighbour/list.hpp"
#include "libatom/potentials/base.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::potentials {

  void Dimer::gradient(SimCell &cell, neighbour::List const &nl, std::size_t num_threads) {
    //

    m_core.clear();

    m_wrapped->gradient(cell, nl, num_threads);  // Gradient at centre (g1)
    m_g0 = cell(Gradient{});

    using std::swap;
    swap(cell(Position{}), m_active);  // Save atoms positions for later.
   
    cell(Position{}) = m_active + m_opt.delta_r * cell(Axis{});

    m_wrapped->gradient(cell, nl, num_threads);  // Gradient at end (g1)
    m_g1 = cell(Gradient{}); 

    floating curv = [&]{
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

                cell(Position{}) = m_active + m_opt.delta_r * m_axisp;

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
                    constexpr auto str= "\tDimer: i={:<4} theta={:f} curv={:f}\n";
                    fmt::print(str, i, theta_min, gdot(m_g1 - m_g0, cell(Axis{})) / m_opt.delta_r);
                }


                if (std::abs(theta_min) < m_opt.theta_tol) {
                    return gdot(m_g1 - m_g0, cell(Axis{})) / m_opt.delta_r;
                }               
            }
        }
    }();

    swap(cell(Position{}), m_active);     // Return cell to original state
                 
    if (m_opt.relax_in_convex && curv > 0) {
        cell(Gradient{}) = -gdot(m_g0, cell(Axis{})) * cell(Axis{});
    } else {
        cell(Gradient{}) = m_g0 - 2 * gdot(m_g0, cell(Axis{})) * cell(Axis{});
    
    }

    // return curv;
  }

}  // namespace otf::potentials