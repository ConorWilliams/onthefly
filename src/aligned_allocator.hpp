#pragma once

#include <cstddef>
#include <new>
#include <type_traits>

// Minimal C++17 and above (over)aligned allocator
template <class T, std::size_t Align> class aligned {
  public:
    // Boilerplate as std::allocator
    using value_type = T;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    using propagate_on_container_move_assignment = std::true_type;
    using is_always_equal = std::true_type;

    // Required due to NTTP
    template <class U> struct rebind { using other = aligned<U, Align>; };

    constexpr aligned() noexcept = default;
    constexpr aligned(aligned const &) noexcept = default;
    template <class U> constexpr aligned(aligned<U, Align> const &) noexcept {}

    [[nodiscard]] static constexpr T *allocate(std::size_t n) {
        return static_cast<T *>(::operator new(n * sizeof(T), align));
    }

    static constexpr void deallocate(T *p, std::size_t) { ::operator delete(p, align); }

  private:
    static_assert(Align > 0 && !(Align & (Align - 1)), "Align must be a power of 2");

    static constexpr std::align_val_t align{Align > alignof(T) ? Align : alignof(T)};

    template <class U>
    friend constexpr bool operator==(aligned const &, aligned<U, Align> const &) noexcept {
        return true;
    }

    template <class U>
    friend constexpr bool operator!=(aligned const &, aligned<U, Align> const &) noexcept {
        return false;
    }
};
