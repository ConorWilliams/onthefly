
#include "catalogue.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cereal/archives/binary.hpp"
#include "cereal/archives/json.hpp"
#include "cereal/archives/portable_binary.hpp"
#include "cereal/archives/xml.hpp"
#include "config.hpp"
#include "local/discrete_key.hpp"
#include "local/environment.hpp"
#include "local/geometry.hpp"
#include "package/package.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

namespace options {

Catalogue Catalogue::load(toml::v2::table const &config) {
    Catalogue opt;

    opt.r_env = fetch<double>(config, "catalogue", "r_env");
    opt.delta = fetch<double>(config, "catalogue", "delta");
    opt.match_best = fetch<bool>(config, "catalogue", "match_best");

    opt.format = config["catalogue"]["format"].value_or(opt.format);
    opt.fname = config["catalogue"]["fname"].value_or(opt.fname);
    opt.load_from_disk = config["catalogue"]["load_from_disk"].value_or(opt.load_from_disk);

    return opt;
}

}  // namespace options

Catalogue::Catalogue(options::Catalogue const &opt) : _opt{opt} {
    //
    if (!_opt.load_from_disk) {
        return;
    }

    std::ifstream file(_opt.fname);

    if (!file.good()) {
        std::cout << "Could not open: \"" << _opt.fname << "\"\n";
        std::cout << "Assuming this is first run, initialising empty catalogue\n";
        return;
    }

    Catalogue cat;

    try {
        if (_opt.format == "binary") {
            cereal::BinaryInputArchive iarchive(file);
            iarchive(cat);
        } else if (_opt.format == "json") {
            cereal::JSONInputArchive iarchive(file);
            iarchive(cat);
        } else if (_opt.format == "portable_binary") {
            cereal::PortableBinaryInputArchive iarchive(file);
            iarchive(cat);
        } else if (_opt.format == "xml") {
            cereal::XMLInputArchive iarchive(file);
            iarchive(cat);
        } else {
            ALWAYS_CHECK(false, "Invalid catalogue.format");
        }
    } catch (...) {
        std::cerr << "Could not load catalogue, mismatched format\n";
        throw;
    }

    ALWAYS_CHECK(_opt.r_env == cat._opt.r_env, "Catalogue incompatible with .r_env");
    ALWAYS_CHECK(_opt.delta == cat._opt.delta, "Catalogue incompatible with .delta");

    _size = std::move(cat._size);
    _catalogue = std::move(cat._catalogue);

    optimise();  // Early call as nothing holds pointers yet
}

void Catalogue::optimise() {
    for (auto &&[k, v] : _catalogue) {
        std::stable_sort(
            v.begin(), v.end(), [](auto const &a, auto const &b) { return a.freq > b.freq; });
    }
}

std::vector<std::size_t> Catalogue::canon_update(std::vector<DiscreteKey> const &keys,
                                                 std::vector<Geometry> &geos,
                                                 std::vector<Catalogue::pointer> &env) {
    env.clear();

    std::vector<std::size_t> out;

    for (std::size_t i = 0; i < keys.size(); ++i) {
        auto &&[ptr, inserted] = canon_try_emplace(keys[i], geos[i]);

        env.push_back(ptr);

        if (ptr->freq++ == 0) {
            out.push_back(i);
        }
    }

    return out;
}

//   Set all frequencies to zero
void Catalogue::reset_counts() {
    for (auto &&[k, v] : _catalogue) {
        for (auto &&env : v) {
            env.freq = 0;
        }
    }
}

void Catalogue::report() const {
    for (auto &&[k, v] : _catalogue) {
        std::cout << k.centre_col << " |";

        for (auto &&c : k.sdf) {
            std::cout << ' ' << c;
        }

        std::cout << " |";

        for (auto &&t : v) {
            std::cout << " " << t.freq;
        }

        std::cout << '\n';
    }

    std::cout << "Found " << _size << " unique topologies in " << _catalogue.size() << " bins!\n";
}

void Catalogue::write() const {
    //

    std::cout << "[[WRITE]]\n";

    std::ofstream file(_opt.fname);

    if (_opt.format == "binary") {
        cereal::BinaryOutputArchive oarchive(file);
        oarchive(*this);
    } else if (_opt.format == "json") {
        cereal::JSONOutputArchive oarchive(file);
        oarchive(*this);
    } else if (_opt.format == "portable_binary") {
        cereal::PortableBinaryOutputArchive oarchive(file);
        oarchive(*this);
    } else if (_opt.format == "xml") {
        cereal::XMLOutputArchive oarchive(file);
        oarchive(*this);
    } else {
        ALWAYS_CHECK(false, "Invalid catalogue.format");
    }
}

template <typename It> It Catalogue::lin_search(It beg, It end, Geometry &mut) const {
    if (_opt.match_best) {
        if (beg == end) {
            return end;
        }

        It first = beg;
        It min = first;

        double min_proj = first->geo.norm(mut);

        for (++first; first != end; ++first) {
            double proj = first->geo.norm(mut);

            if (proj < min_proj) {
                min = first;
                min_proj = proj;
            }
        }

        if (mut.permute_onto(min->delta, min->geo)) {
            return min;
        } else {
            return end;
        }
    } else {
        return std::find_if(beg, end, [&](Environment const &ref) -> bool {
            // Test if fuzzy keys match (fast)
            if (!ref.geo.equiv(ref.delta, mut)) {
                return false;
            }

            // Full-Monte equivalence
            return static_cast<bool>(mut.permute_onto(ref.delta, ref.geo));
        });
    }
}

// Converts geo into canonical order, inserts into _catalogue if not already there and returns a
// reference to the topology equivalent to "t" in _catalogue.
std::pair<Catalogue::pointer, bool> Catalogue::canon_try_emplace(DiscreteKey const &key,
                                                                 Geometry &geo) {
    // "it" always points to valid bucket (possibly empty)
    auto [it, inserted] = _catalogue.try_emplace(key);

    if (!inserted) {
        // Existing key, must search bucket for explicit match;
        auto match = lin_search(it->second.begin(), it->second.end(), geo);

        // If found a match, return it
        if (match != it->second.end()) {
            return {pointer(it, match - it->second.begin()), false};
        }
    }

    // Otherwise insert new geo at end of bucket
    it->second.emplace_back(geo, _opt.delta);
    ++_size;

    return {pointer(it, it->second.size() - 1), true};
}

bool Catalogue::try_canon(std::vector<DiscreteKey> const &keys,
                          std::vector<Geometry> &geos,
                          std::vector<Catalogue::pointer> &env) {
    env.clear();

    for (std::size_t i = 0; i < keys.size(); ++i) {
        //
        auto [it, inserted] = _catalogue.try_emplace(keys[i]);

        if (!inserted) {
            // Existing key, must search bucket for explicit match;
            auto match = lin_search(it->second.begin(), it->second.end(), geos[i]);

            // If found a match, return it
            if (match != it->second.end()) {
                env.push_back(pointer{it, match - it->second.begin()});
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    return true;
}
