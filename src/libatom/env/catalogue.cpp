
#include "libatom/env/catalogue.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/env/geometry.hpp"
#include "libatom/env/topology.hpp"
#include "libatom/utils.hpp"

namespace otf::env {

  std::optional<Catalogue::Pointer> Catalogue::canon_find(LocalEnv& mut) {
    // "it" always points to valid bucket (possibly empty)
    auto [it, inserted] = m_cat.try_emplace(mut.key());

    if (!inserted) {
      // Existing key, must search bucket for explicit match
      auto match = std::find_if(it->second.begin(), it->second.end(), [&](Environment const& ref) {
        //
        floating r_min = std::min(mut.fingerprint().r_min(), ref.fingerprint.r_min());

        floating delta = std::min(0.4 * ref.delta_mod * r_min, m_opt.delta_max);

        // fmt::print(stderr, "Delta:{}, r_min:{}\n", delta, r_min);

        // Test if fuzzy keys match (fast)
        if (!equiv(ref.fingerprint, mut.fingerprint(), M_SQRT2 * delta)) {
          return false;
        }

        if (mut.permute_onto(ref.ref_geo, delta)) {
          return true;
        } else {
          //   fmt::print(stderr, "False equiv\n");
          return false;
        }
      });

      // If found a match, return it
      if (match != it->second.end()) {
        return Pointer(it, match - it->second.begin());
      }
    }

    return std::nullopt;
  }

  Catalogue::Pointer Catalogue::insert(LocalEnv const& env) {
    //
    // ASSERT(!canon_find(env), "Environemnt already in catalogue.");

    auto [it, inserted] = m_cat.try_emplace(env.key());

    auto& new_ref = it->second.emplace_back(env.fingerprint());

    for (auto&& elem : env) {
      new_ref.ref_geo.emplace_back(elem(Position{}), elem(Colour{}));
    }

    m_size++;

    return {it, it->second.size() - 1};
  }

}  // namespace otf::env