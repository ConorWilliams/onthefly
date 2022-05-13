#pragma once

#include <fstream>  // Required by toml++

#include "indicators/block_progress_bar.hpp"
#include "indicators/cursor_control.hpp"
#include "indicators/progress_bar.hpp"
#include "kinetics/basin.hpp"
#include "kinetics/superbasin.hpp"
#include "kinetics/supercache.hpp"
#include "local/catalogue.hpp"
#include "local/classify.hpp"
#include "local/geometry.hpp"
#include "minimise/minimiser_base.hpp"
#include "package/package.hpp"
#include "potentials/potential_base.hpp"
#include "sp_search/find_mech.hpp"
#include "sp_search/sp_search_base.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"
#include "visualise.hpp"

// RAII wrapper
struct Bar : indicators::ProgressBar {
    Bar(std::size_t size)
        : indicators::ProgressBar{indicators::option::BarWidth{50},
                                  indicators::option::Start{" ["},
                                  indicators::option::Fill{"="},
                                  indicators::option::Lead{">"},
                                  indicators::option::Remainder{"-"},
                                  indicators::option::End{"]"},
                                  indicators::option::PrefixText{std::to_string(size) + " iters:"},
                                  indicators::option::ShowElapsedTime{true},
                                  indicators::option::ShowRemainingTime{true},
                                  indicators::option::MaxProgress{size}} {
        indicators::show_console_cursor(false);
        set_progress(0);
    }

    ~Bar() { indicators::show_console_cursor(true); }
};
