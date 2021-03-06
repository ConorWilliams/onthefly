
#pragma once

#include <string>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/utils.hpp"

namespace otf::potentials {

/**
 * @brief Computes a set of natural cubic spline coefficients for a uniformly tabulated function.
 *
 * See Wikipedia algorithm: https://en.wikipedia.org/wiki/Spline_(mathematics)
 */
class Spline {
public:
  /**
   * @brief Construct an empty Spline
   */
  Spline() = default;

  /**
   * @brief Construct a new Spline object from (n + 1) evenly spaced tabulated values on the
   * interval {0, dx, ..., ndx}.
   */
  Spline(std::vector<floating> y, floating dx);

  /**
   * @brief Interpolate tabulated function.
   *
   * @param x
   * @return f(x)
   */
  [[nodiscard]] floating f(floating x) const {
    auto [x0, spine] = fetch(x);
    return spine.a + x0 * (spine.b + x0 * (spine.c + x0 * spine.d));
  }

  /**
   * @brief Interpolate tabulated functions gradient.
   *
   * @param x
   * @return f'(x)
   */
  [[nodiscard]] floating fp(floating x) const {
    auto [x0, spine] = fetch(x);
    return spine.b + x0 * (2.0 * spine.c + x0 * 3.0 * spine.d);
  }

  /**
   * @brief Interpolate tabulated second derivative.
   *
   * Note this may not be continuous or smooth.
   *
   * @param x
   * @return f''(x)
   */
  [[nodiscard]] floating fpp(floating x) const {
    auto [x0, spine] = fetch(x);
    return 2 * spine.c + 6.0 * x0 * spine.d;
  }

private:
  struct Spine {
    floating a, b, c, d;
  };

  std::vector<Spine> m_spines;

  floating m_dx = 0;
  floating m_inv_dx = 0;

  std::pair<floating, Spine const&> fetch(floating x) const {
    ASSERT(x >= 0, "X is less than zero");
    std::size_t i = static_cast<std::size_t>(x * m_inv_dx);
    ASSERT(i < m_spines.size(), "X is outside tabulated region");
    return {x - i * m_dx, m_spines[i]};
  }
};

} // namespace otf::potentials