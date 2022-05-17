#include "libatom/utils.hpp"

#include <string_view>

#include "doctest/doctest.h"
#include "onthefly/version.hpp"

TEST_CASE("OnTheFly version") {
  static_assert(std::string_view(ONTHEFLY_VERSION) == std::string_view("0.1"));
}

TEST_CASE("ipow") {
  //
  CHECK(otf::ipow<0>(10) == 1);
  CHECK(otf::ipow<1>(10) == 10);
  CHECK(otf::ipow<2>(10) == 100);
  CHECK(otf::ipow<3>(10) == 1000);
  CHECK(otf::ipow<4>(10) == 10000);
  CHECK(otf::ipow<5>(10) == 100000);
}
