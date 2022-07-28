#pragma once

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

struct Mechanism {};

class Catalogue {
public:
  struct Options {
    /** @brief Maximum difference in L2 norm between LEs for them to be considered the same. */
    double delta_max = 0.5;
    /** @brief Smaller to decrease false positives but, values < 1.0 introduce false negatives. */
    double overfuzz = 0.5;
  };

  /**
   * @brief Struct used to represent LocalEnv in the catalogue
   *
   */
  struct Env : Geometry<Position, Colour> {
    //
    LocalEnv::Fingerprint fingerprint;

    std::vector<Mechanism> mechs{};
    std::size_t freq = 1;
    floating delta_mod = 1;

    Env() = default;

    explicit Env(LocalEnv::Fingerprint const& f) : fingerprint(f) {}
  };

  /**
   * @brief  Only invalidated if ...
   */
  struct Pointer {
  public:
    Env* operator->() const noexcept {
      return m_it->second.data() + m_offset;
    }

  private:
    friend class Catalogue;

    Pointer() = default;

    using iterator = std::map<LocalEnv::Key, std::vector<Env>, std::less<>>::iterator;

    Pointer(iterator it, std::size_t offset) noexcept : m_it{it}, m_offset{offset} {}

    iterator m_it;
    std::size_t m_offset;
  };

  explicit Catalogue(Options const& opt) : m_opt(opt) {}

  std::size_t size() const noexcept {
    return m_size;
  }

  std::size_t num_keys() const noexcept {
    return m_cat.size();
  }

  void optimize() {
    for (auto&& [k, v] : m_cat) {
      std::sort(v.begin(), v.end(), [](Env const& e1, Env const& e2) {
        // Biggest first
        return e1.freq > e2.freq;
      });
    }
  }

  /**
   * @brief Simultaniously canonise the input local environment and find its match in the
   * catalogue.
   *
   * This is a deterministic search and will find the first match that is equivilent to env.
   *
   * Each call to canon_find increases the frequency count if a match is found.
   *
   * @return std::optional<Pointer> Empty if no match.
   */
  [[nodiscard]] std::optional<Pointer> canon_find(LocalEnv& env);

  /**
   * @brief Simultaniously canonise the input local environment and insert it into the Catalogue.
   *
   * Should only be called if env is not already in the catalogue.
   */
  Pointer insert(LocalEnv const& env);

private:
  Options m_opt;
  std::size_t m_size = 0;
  std::map<LocalEnv::Key, std::vector<Env>, std::less<>> m_cat;
};

} // namespace otf::env