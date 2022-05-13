#include "sp_search/sp_search_base.hpp"

#include <memory>
#include <stdexcept>
#include <string>

#include "config.hpp"
#include "sp_search/dimer/dimer.hpp"
#include "sp_search/dimer/l_shrink.hpp"
#include "sp_search/dimer/search.hpp"
#include "sp_search/dimer/shrinking.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

// Customisation point, dynamically select minimiser
std::unique_ptr<SearchBase> load_sp_search(toml::v2::table const& config) {
    //
    std::string kind = fetch<std::string>(config, "sp_search", "kind");

    if (kind == "Dimer") {
        return std::make_unique<DimerSPS<Dimer>>(load_minimiser(config),
                                                 Dimer{options::Dimer::load(config)});
    } else if (kind == "LShrink") {
        return std::make_unique<DimerSPS<LShrink>>(load_minimiser(config),
                                                   LShrink{options::LShrink::load(config)});
    } else if (kind == "Shrinking") {
        return std::make_unique<DimerSPS<ShrinkingDimer>>(
            load_minimiser(config), ShrinkingDimer{options::ShrinkingDimer::load(config)});
    } else {
        throw std::runtime_error("Unsupported minimise field selected : " + kind);
    }
}