#pragma once

#include <cstdint>

#include "libatom/asserts.hpp"
#include "libatom/system/atomvector.hpp"
#include "libatom/system/simbox.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief A SimCell is a collection of active/frozen atoms and inherits from an OrthoSimBox
   *
   */
  class SimCell : public OrthoSimBox {
  private:
    AtomVector m_active;
    AtomVector m_frozen;

  public:
    explicit SimCell(OrthoSimBox const& box) : OrthoSimBox{box} {}

    /**
     * @brief Get the total (active + frozen) number of atoms in the SimCell
     */
    std::size_t size() const noexcept { return m_active.size() + m_frozen.size(); }

    /**
     * @brief Get a const reference to the active atoms
     */
    AtomVector const& active() const noexcept { return m_active; }

    /**
     * @brief Get a const reference to the frozen atoms
     */
    AtomVector const& frozen() const noexcept { return m_frozen; }

    /**
     * @brief Get a reference to the active atoms
     */
    AtomVector& active() noexcept { return m_active; }

    /**
     * @brief Get a reference to the frozen atoms
     */
    AtomVector& frozen() noexcept { return m_frozen; }

    // VecN<double> active_disp(VecN<double> const &others) const;
  };

}  // namespace otf