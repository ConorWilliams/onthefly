
#include "potentials/potential_base.hpp"

#include <memory>
#include <stdexcept>
#include <string>

#include "config.hpp"
#include "potentials/ADP/potential.hpp"
#include "potentials/EAM/potential.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

// Customisation point, dynamically select potentials
std::unique_ptr<PotentialBase> load_potential(toml::v2::table const& config) {
    //
    std::string kind = fetch<std::string>(config, "potential", "kind");

    if (kind == "EAM") {
        return std::make_unique<PotentialEAM>(load_eam(config));
    } else if (kind == "ADP") {
        return std::make_unique<PotentialADP>(load_adp(config));
    } else {
        throw std::runtime_error("Unsupported potential selected : " + kind);
    }
}
