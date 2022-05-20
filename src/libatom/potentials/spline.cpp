#include "libatom/potentials/spline.hpp"

#include <cstddef>

#include "libatom/utils.hpp"

namespace otf {

  Spline::Spline(std::vector<floating> y, floating dx) : m_spines{}, m_dx(dx), m_inv_dx(1 / dx) {
    //
    // Extend spline with one value such that if called with x = ndx which rounds up we have no
    // issues.
    y.push_back(floating(y.back()));

    std::size_t n = y.size() - 1;

    // 1
    std::vector<floating> a(n + 1);
    // 2
    std::vector<floating> b(n);
    std::vector<floating> d(n);
    // 4
    std::vector<floating> alpha(n);
    // 5
    std::vector<floating> c(n + 1);
    std::vector<floating> l(n + 1);
    std::vector<floating> mu(n + 1);
    std::vector<floating> z(n + 1);

    // 1
    for (std::size_t i = 0; i <= n; ++i) {
      a[i] = y[i];
    }

    // 3
    // h_i = _dx

    // 4
    for (std::size_t i = 1; i <= n - 1; ++i) {
      const floating inv_h = 3 / dx;

      alpha[i] = inv_h * (a[i + 1] - 2 * a[i] + a[i - 1]);
    }

    // 6
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    // 7
    for (std::size_t i = 1; i <= n - 1; ++i) {
      l[i] = dx * (4 - mu[i - 1]);
      mu[i] = dx / l[i];
      z[i] = (alpha[i] - dx * z[i - 1]) / l[i];
    }

    // 8
    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    // 9
    for (std::size_t i = 1; i <= n; ++i) {
      std::size_t j = n - i;

      c[j] = z[j] - mu[j] * c[j + 1];
      b[j] = (a[j + 1] - a[j]) / dx - dx * (c[j + 1] + 2 * c[j]) / 3;
      d[j] = (c[j + 1] - c[j]) / (3 * dx);
    }

    // 11
    for (std::size_t i = 0; i <= n - 1; ++i) {
      m_spines.push_back({a[i], b[i], c[i], d[i]});
    }
  }

}  // namespace otf