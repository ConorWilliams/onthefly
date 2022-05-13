

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <type_traits>

#include "doctest/doctest.h"
#include "fmt/core.h"
#include "libatom/asserts.hpp"
#include "libatom/utils.hpp"

template <typename T, std::size_t Extent = 1> struct Member {
  //
  using scalar_type = T;
  using vector_type = std::conditional_t<Extent == 1, T, Eigen::Array<T, Extent, 1>>;

  static constexpr std::size_t extent = Extent;
};

template <typename Tag> class EigenArrayAdaptor {
private:
  static_assert(std::is_empty_v<Tag>, "Tag types are required to be empty");

  Eigen::Array<typename Tag::scalar_type, Tag::extent, Eigen::Dynamic> m_data;

public:
  auto const &raw_array(Tag) const { return m_data; };

  auto &raw_array(Tag) { return m_data; };

  decltype(auto) get(Tag, std::size_t i) const {
    if constexpr (Tag::extent == 1) {
      return m_data[i];
    } else {
      return m_data.col(i);
    }
  };

  decltype(auto) get(Tag, std::size_t i) {
    if constexpr (Tag::extent == 1) {
      return m_data[i];
    } else {
      return m_data.col(i);
    }
  };
};

template <typename... Mems> class AtomArray : private EigenArrayAdaptor<Mems>... {
private:
  std::size_t m_size = 0;
  std::size_t m_cap = 0;

  using EigenArrayAdaptor<Mems>::raw_array...;
  using EigenArrayAdaptor<Mems>::get...;

public:
  std::size_t size() const { return size; }

  /**
   * @brief Destructive resize, destroyes all atom data, allocates unless size() == size.
   */
  void resize(std::size_t size) {
    m_size = size;
    m_cap = size;

    ((raw_array(Mems{}).resize(Eigen::NoChange, m_cap)), ...);
  }

  template <typename Tag> decltype(auto) operator()(Tag, std::size_t i) const {
    STACK();
    ASSERT(i < m_size, "Out of bounds");
    return get(Tag{}, i);
  }

  template <typename Tag> decltype(auto) operator()(Tag, std::size_t i) {
    STACK();
    ASSERT(i < m_size, "Out of bounds");
    return get(Tag{}, i);
  }

  template <typename Tag> auto operator()(Tag) const { return raw_array(Tag{}).leftCols(m_size); }

  template <typename Tag> auto operator()(Tag) { return raw_array(Tag{}).leftCols(m_size); }

  void push_back(typename Mems::vector_type const &...args) { emplace_back(args...); }

  template <typename... Args> void emplace_back(Args &&...args) {
    //
    static_assert(sizeof...(Args) == sizeof...(Mems), "Too few args to emplace_back");

    STACK();

    ASSERT(m_size <= m_cap, "Pre check invariant");

    if (m_size == m_cap) {
      // Need to grow arrays
      m_cap = std::max(m_cap + 1, std::size_t(1.5 * m_cap));

      ((raw_array(Mems{}).conservativeResize(Eigen::NoChange, m_cap)), ...);
    }

    ((get(Mems{}, m_size) = std::forward<Args>(args)), ...);

    m_size += 1;

    ASSERT(m_size <= m_cap, "Post check invariant");
  }
};

TEST_CASE("avec") {
  //
  struct Pos : Member<double, 3> {};
  struct State : Member<int> {};
  struct AtNum : Member<std::vector<int>> {};

  AtomArray<Pos, State, AtNum> arr{};

  arr.push_back({1, 1, 1}, 2, {1, 2, 3});

  for (auto &&elem : arr(AtNum{}, 0)) {
    fmt::print("{}", elem);
  }
}