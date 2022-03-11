#include "onthefly/onthefly.hpp"

#include <string>

#include "doctest/doctest.h"
#include "onthefly/version.h"

TEST_CASE("OnTheFly") {
  using namespace onthefly;

  OnTheFly onthefly("Tests");

  CHECK(onthefly.greet(LanguageCode::EN) == "Hello, Tests!");
  CHECK(onthefly.greet(LanguageCode::DE) == "Hallo Tests!");
  CHECK(onthefly.greet(LanguageCode::ES) == "¡Hola Tests!");
  CHECK(onthefly.greet(LanguageCode::FR) == "Bonjour Tests!");
}

TEST_CASE("OnTheFly version") {
  static_assert(std::string_view(ONTHEFLY_VERSION) == std::string_view("1.0"));
  CHECK(std::string(ONTHEFLY_VERSION) == std::string("1.0"));
}
