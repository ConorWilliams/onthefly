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
    [[nodiscard]] Position::matrix_type& newton_step(SimCell const&);

  private:
    std::size_t m_n;
    std::size_t m_k;

    struct Elem {
      Position::matrix_type s;
      Gradient::matrix_type y;
      floating rho;
      floating alpha;
    };

    std::vector<Elem> m_hist;

    Position::matrix_type m_prev_x;
    Gradient::matrix_type m_prev_g;

    Position::matrix_type m_r;
    Gradient::matrix_type m_q;
  };

}  // namespace otf::minimise