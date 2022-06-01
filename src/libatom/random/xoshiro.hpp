#pragma once

// The code in this file is adapted from the original implementation:
// http://prng.di.unimi.it/xoshiro256starstar.c

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "libatom/asserts.hpp"

namespace otf::random {

  /**
   * @brief Xoshiro utility frunction.
   */
  [[nodiscard]] inline constexpr uint64_t rotl(uint64_t const x, int const k) noexcept {
    return (x << k) | (x >> (64 - k));
  }

  /**
   * @brief This is a <random> compatable implementation of the xoshiro256** 1.0 PRNG
   *
   * From the origional: "An all-purpose, rock-solid generator. It has excellent (sub-ns) speed, a
   * state (256 bits) that is large enough for any parallel application, and it passes all tests we
   * are aware of. The state must be seeded so that it is not everywhere zero."
   */
  class Xoshiro {
  public:
    explicit constexpr Xoshiro(std::array<std::uint64_t, 4> state) noexcept : m_state{state} {
      for ([[maybe_unused]] auto&& elem : m_state) {
        ASSERT(elem != 0, "Bad seed!");
      }
    }

    /**
     * @brief Get the minimum value of the generator.
     */
    [[nodiscard]] static constexpr std::uint64_t min() noexcept {
      return std::numeric_limits<std::uint64_t>::lowest();
    }

    /**
     * @brief Get the maximum value of the generator.
     */
    [[nodiscard]] static constexpr std::uint64_t max() noexcept {
      return std::numeric_limits<std::uint64_t>::max();
    }

    /**
     * @brief Generate a random bit sequence and advance the state of the generator.
     */
    constexpr std::uint64_t operator()() noexcept {
      std::uint64_t const result = rotl(m_state[1] * 5, 7) * 9;

      std::uint64_t const t = m_state[1] << 17;

      m_state[2] ^= m_state[0];
      m_state[3] ^= m_state[1];
      m_state[1] ^= m_state[2];
      m_state[0] ^= m_state[3];

      m_state[2] ^= t;

      m_state[3] = rotl(m_state[3], 45);

      return result;
    }

    /**
     * @brief This is the jump function for the generator.
     *
     * It is equivalent to 2^128 calls to operator(); it can be used to generate 2^128
     * non-overlapping subsequences for parallel computations.
     */
    constexpr void jump() noexcept {
      constexpr std::uint64_t JUMP[] = {
          0x180ec6d33cfd0aba,
          0xd5a61266f0c9392c,
          0xa9582618e03fc9aa,
          0x39abdc4529b1661c,
      };

      std::uint64_t s0 = 0;
      std::uint64_t s1 = 0;
      std::uint64_t s2 = 0;
      std::uint64_t s3 = 0;

      for (std::uint64_t jump : JUMP) {
        for (int b = 0; b < 64; ++b) {
          if (jump & UINT64_C(1) << b) {
            s0 ^= m_state[0];
            s1 ^= m_state[1];
            s2 ^= m_state[2];
            s3 ^= m_state[3];
          }
          operator()();
        }
      }

      m_state[0] = s0;
      m_state[1] = s1;
      m_state[2] = s2;
      m_state[3] = s3;
    }

    /**
     * @brief This is the long-jump function for the generator.
     *
     * It is equivalent to 2^192 calls to operator(); it can be used to generate 2^64 starting
     * points, from each of which jump() will generate 2^64 non-overlapping subsequences for
     * parallel distributed computations.
     */
    constexpr void long_jump() noexcept {
      constexpr std::uint64_t LONG_JUMP[] = {
          0x76e15d3efefdcbbf,
          0xc5004e441c522fb3,
          0x77710069854ee241,
          0x39109bb02acbe635,
      };

      std::uint64_t s0 = 0;
      std::uint64_t s1 = 0;
      std::uint64_t s2 = 0;
      std::uint64_t s3 = 0;

      for (std::uint64_t jump : LONG_JUMP) {
        for (int b = 0; b < 64; ++b) {
          if (jump & UINT64_C(1) << b) {
            s0 ^= m_state[0];
            s1 ^= m_state[1];
            s2 ^= m_state[2];
            s3 ^= m_state[3];
          }
          operator()();
        }
      }
      m_state[0] = s0;
      m_state[1] = s1;
      m_state[2] = s2;
      m_state[3] = s3;
    }

  private:
    std::array<std::uint64_t, 4> m_state;
  };

}  // namespace otf::random
