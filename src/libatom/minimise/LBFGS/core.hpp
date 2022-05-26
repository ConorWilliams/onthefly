#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/atom_array.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::minimise {

  /**
   * @brief Holds variable history and computes the next lbfgs newton step at each call.
   */
  class CoreLBFGS {
  public:
    CoreLBFGS(std::size_t n) : m_n{n}, m_k{0}, m_hist(n){};

    /**
     * @brief Reset history for re-use.
     */
    void clear() { m_k = 0; }

    /**
     * @brief Computes the product of the approximate hessian and gradiant, H ∇f, using the L-BFGS
     * two-loop recursion.
     *
     * The newton step towards the minimum is x -= a * H ∇f with a = 1 the best guess.
     *
     * @return A view of the newton step array, H ∇f, (it is ok to modify this view it will be
     * overwritten upon next call).
     */
    [[nodiscard]] SimCell::underlying_t<Gradient>& newton_step(SimCell const&);

  private:
    std::size_t m_n;
    std::size_t m_k;

    struct Elem {
      SimCell::underlying_t<Position> s;
      SimCell::underlying_t<Gradient> y;
      floating rho;
      floating alpha;
    };

    std::vector<Elem> m_hist;

    SimCell::underlying_t<Position> m_prev_x;
    SimCell::underlying_t<Gradient> m_prev_g;

    SimCell::underlying_t<Gradient> m_q;
    SimCell::underlying_t<Gradient> m_r;
  };

}  // namespace otf::minimise