#include "minimise/minimiser_base.hpp"

#include <memory>
#include <stdexcept>
#include <string>

#include "config.hpp"
#include "minimise/BB/bb.hpp"
#include "minimise/HYBRID/hybrid.hpp"
#include "minimise/LBFGS/lbfgs.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

// Customisation point, dynamically select minimiser
std::unique_ptr<MinimiserBase> load_minimiser(toml::v2::table const& config) {
    //
    std::string kind = fetch<std::string>(config, "minimiser", "kind");

    if (kind == "BB") {
        return std::make_unique<MinimiseBB>(options::MinimiseBB::load(config));
    } else if (kind == "HYBRID") {
        return std::make_unique<MinimiseHYBRID>(options::MinimiseHYBRID::load(config));
    } else if (kind == "LBFGS") {
        return std::make_unique<MinimiseLBFGS>(options::MinimiseLBFGS::load(config));
    } else {
        throw std::runtime_error("Unsupported minimise field selected : " + kind);
    }
}