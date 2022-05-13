#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

#include "cereal/types/array.hpp"
#include "config.hpp"

// Error handler, no inline for smaller binaries
[[noreturn]] __attribute__((noinline)) void assert_handler(int lineno,
                                                           char const *file,
                                                           char const *condition,
                                                           std::string message);

// Custom assert-like macros with error messages
#define ALWAYS_CHECK(condition, message)                             \
    do {                                                             \
        if (!(condition)) {                                          \
            assert_handler(__LINE__, __FILE__, #condition, message); \
        }                                                            \
    } while (false)

#if defined(NDEBUG) || defined(OLKMC_NO_DEBUG)
#    define CHECK(condition, message) \
        do {                          \
        } while (false)
#else
#    define CHECK(condition, message) ALWAYS_CHECK(condition, message)
#endif

// Throwing version of std::getline
std::istringstream safe_getline(std::ifstream &file);

template <typename T, typename Tbl> std::optional<T> try_fetch(Tbl const &tbl) {
    return tbl.template value_exact<T>();
}

template <typename T, typename Tbl, typename... Args>
std::optional<T> try_fetch(Tbl const &tbl, std::string_view str, Args &&...args) {
    return try_fetch<T>(tbl[str], std::forward<Args>(args)...);
}

// Helper to extract values from config while throwing usefull errors
template <typename T, typename Tbl, typename... Args>
T fetch(Tbl const &tbl, std::string_view str, Args &&...args) {
    //
    if (std::optional<T> maybe = try_fetch<T>(tbl[str], std::forward<Args>(args)...)) {
        return *maybe;
    } else {
        std::stringstream ss;

        ss << "Config missing: " << str;

        ((ss << '.' << std::forward<Args>(args)), ...);

        ALWAYS_CHECK(false, ss.str());
    }
}

// Timing utilities for quick testing

struct clock_tick {
    std::string name;
    typename std::chrono::high_resolution_clock::time_point start;
};

inline clock_tick tick(std::string const &name, bool print = false) {
    if (print) {
        std::cout << "Timing: " << name << '\n';
    }
    return {name, std::chrono::high_resolution_clock::now()};
}

template <typename... Args> int tock(clock_tick &x, Args &&...args) {
    using namespace std::chrono;

    auto const stop = high_resolution_clock::now();

    auto const time = duration_cast<microseconds>(stop - x.start).count();

    std::cout << x.name << ": " << time << "/us";

    (static_cast<void>(std::cout << ',' << ' ' << args), ...);

    std::cout << std::endl;

    return time;
}

// L2 squared-norm, generic eigen-array support
template <typename E> double norm_sq(E const &dr) { return (dr * dr).sum(); }

// L2 norm, generic eigen-array support
template <typename E> double norm(E &&dr) { return std::sqrt(norm_sq(std::forward<E>(dr))); }

// Dot product for generic eigen-arrays
template <typename A, typename B> double dot(A &&a, B &&b) {
    return (std::forward<A>(a) * std::forward<B>(b)).sum();
}

// Specialise Eigen types for cereal
namespace cereal {

template <class Archive, class T> void serialize(Archive &archive, Vec3<T> &v) {
    archive(v[0], v[1], v[2]);
}

template <class Archive, class T> void serialize(Archive &archive, VecN<T> &v) {
    for (auto &&elem : v) {
        archive(elem);
    }
}

}  // namespace cereal

// Small compile time N by N matrix such that: SymMat a; a(i,j) == a(j,i)
template <typename T, std::size_t N> class SymMat {
  public:
    constexpr T const &operator()(std::size_t i, std::size_t j) const {
        return _data[symmetrise(i, j)];
    }
    constexpr T &operator()(std::size_t i, std::size_t j) { return _data[symmetrise(i, j)]; }

    template <typename F> void foreach (F &&f) {
        for (auto &&elem : _data) {
            f(elem);
        }
    }

    template <class Archive> void serialize(Archive &ar) { ar(_data); }

  private:
    std::array<T, N *(N + 1) / 2> _data;

    static constexpr std::size_t symmetrise(std::size_t i, std::size_t j) {
        CHECK(i < N && j < N, "Bad access");
        return std::max(i, j) + (2 * N - 1 - std::min(i, j)) * std::min(i, j) / 2;
    }
};
