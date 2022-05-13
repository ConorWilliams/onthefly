#pragma once

#include <memory>
#include <optional>

#include "config.hpp"
#include "potentials/potential_base.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"

// Here we define the virtual-interface for saddle-point searchers in OLKMC,
class SearchBase {
  public:
    // Make a copy of this
    virtual std::unique_ptr<SearchBase> clone() const = 0;

    // `init` is the initial basin of the potential, `dimer` is the start point for the SP-search
    // final will be the adjacent minima if the search is successful.
    virtual bool find_sp(Supercell const& init,
                         Supercell& dimer,
                         Supercell& final,
                         std::unique_ptr<PotentialBase>& ff)
        = 0;

    // Call parent destructor
    virtual ~SearchBase() {}

  protected:
    // Protected constructor as this is an interface class
    constexpr SearchBase() noexcept = default;
};

// Customisation point, dynamically select minimiser
std::unique_ptr<SearchBase> load_sp_search(toml::v2::table const& config);