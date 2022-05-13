#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "cereal/types/map.hpp"
#include "cereal/types/string.hpp"
#include "config.hpp"
#include "local/discrete_key.hpp"
#include "local/environment.hpp"
#include "local/geometry.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

namespace options {

struct Catalogue {
    double r_env;     // (Angstrom), Radius of local environment
    double delta;     // (Angstrom), Maximum difference in L2 norm between local-environments
    bool match_best;  // If true then selects best instead of first match in catalogue

    std::string format = "portable_binary";
    std::string fname = "olkmc.cat";

    bool load_from_disk = false;

    static Catalogue load(toml::v2::table const &config);

    template <class Archive> void serialize(Archive &ar) { ar(r_env, delta, format, fname); }
};

}  // namespace options

// Data structure to hold the mechanism catalogue
class Catalogue {
  public:
    struct pointer {
        // Invalidated if _catalogue's iterators are invalidated or a call to .optimise() is made
      public:
        Environment *operator->() const {
            CHECK(_bucket_offset >= 0, "Negative offset");
            CHECK(std::size_t(_bucket_offset) < _it->second.size(), "Invalid Offset");
            return _it->second.data() + _bucket_offset;
        }

        // Shrink the environment such that it no longer matches geo and insert the geo as a new
        // environment
        void refine(Geometry &geo) const {
            //
            std::optional<Geometry::Result> perm = geo.permute_onto((*this)->delta, (*this)->geo);

            ALWAYS_CHECK(perm, "refining unmatched");

            std::cout << "Distance of problem child is" << perm->dr << " vs " << (*this)->delta;

            (*this)->delta = perm->dr / 2;

            _it->second.emplace_back(geo, (*this)->delta);
        }

      private:
        friend class Catalogue;

        using map_t = std::map<DiscreteKey, std::vector<Environment>>;

        pointer() = default;

        pointer(map_t::iterator it, std::ptrdiff_t offset) : _it{it}, _bucket_offset{offset} {}

        map_t::iterator _it{};
        std::ptrdiff_t _bucket_offset{};
    };

    Catalogue() = default;

    Catalogue(options::Catalogue const &opt);

    //   Sort buckets into descending frequency order
    void optimise();

    // Update catalogue with new local-environments and sort geos into canonical order, returns
    // indices of atoms in new environments.
    std::vector<std::size_t> canon_update(std::vector<DiscreteKey> const &keys,
                                          std::vector<Geometry> &geos,
                                          std::vector<pointer> &env);

    // Cannonise all geos, if new geo operation fails (returns false)
    bool try_canon(std::vector<DiscreteKey> const &keys,
                   std::vector<Geometry> &geos,
                   std::vector<Catalogue::pointer> &env);

    void write() const;  // To disk

    void reset_counts();

    void report() const;

    template <class Archive> void serialize(Archive &ar) { ar(_opt, _size, _catalogue); }

  private:
    options::Catalogue _opt{};

    std::size_t _size{};  // Number of LEs
    std::map<DiscreteKey, std::vector<Environment>> _catalogue{};

    //////////////////////////////////////////////////////////////////////////////////////////

    // Find Topology in "bucket" that is equivalent to "mut", if match is found then "mut" is
    // permuted on to it
    template <typename It> It lin_search(It beg, It end, Geometry &mut) const;

    // Converts geo into canonical order, inserts into _catalogue if not already there and
    // returns a reference to the topology equivalent to "t" in _catalogue.
    std::pair<pointer, bool> canon_try_emplace(DiscreteKey const &key, Geometry &geo);
};
