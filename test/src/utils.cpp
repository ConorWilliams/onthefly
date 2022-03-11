#include "libatom/utils.hpp"

#include <string_view>

#include "doctest/doctest.h"
#include "onthefly/version.hpp"

TEST_CASE("OnTheFly version") {
  static_assert(std::string_view(ONTHEFLY_VERSION) == std::string_view("0.1"));
}

TEST_CASE("Utils constants version") {
  using namespace otf;

  static_assert(otf::spatial_dims == 2 || otf::spatial_dims == 3);
}
