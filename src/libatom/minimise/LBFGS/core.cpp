
#include "libatom/minimise/LBFGS/core.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::minimise {

  Position::matrix_type& CoreLBFGS::newton_step(SimCell const& atoms) {
    //
    std::size_t prev = (m_k - 1) % m_n;

    // Compute the k-1 th y, s and rho
    if (m_k > 0) {
      //
      m_hist[prev].s = atoms(Position{}) - m_prev_x;
      m_hist[prev].y = atoms(Gradient{}) - m_prev_g;

      // If Wolfie conditions fulfiled during the line search then dot(y, s) > 0. Otherwise we take
      // absolute value to prevent ascent direction.
      m_hist[prev].rho = 1.0 / std::abs(gdot(m_hist[prev].s, m_hist[prev].y));
    }

    m_prev_x = atoms(Position{});
    m_prev_g = atoms(Gradient{});

    m_q = atoms(Gradient{});

    int incur = m_k <= m_n ? 0 : m_k - m_n;
    int bound = m_k <= m_n ? m_k : m_n;

    // Loop 1: over most recent steps
    for (int i = bound - 1; i >= 0; --i) {
      int j = (i + incur) % m_n;

      m_hist[j].alpha = m_hist[j].rho * gdot(m_hist[j].s, m_q);
      m_q -= m_hist[j].alpha * m_hist[j].y;
    }

    // Scaling Hessian_0.
    if (m_k > 0) {
      m_r = m_q * (1.0 / m_hist[prev].rho / gdot(m_hist[prev].y, m_hist[prev].y));
    } else {
      m_r = m_q;  // Start with identity hessian.
    }

    // Loop 2:
    for (int i = 0; i <= bound - 1; ++i) {
      int j = (i + incur) % m_n;

      floating beta = m_hist[j].rho * gdot(m_hist[j].y, m_r);

      m_r += (m_hist[j].alpha - beta) * m_hist[j].s;
    }

    ++m_k;

    return m_r;
  }

}  // namespace otf::minimise