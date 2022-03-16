#include <iostream>

#include "libatom/asserts.hpp"
#include "libatom/external/current_function.hpp"

int test() {
  STACK();

  ASSERT(1 + 1 > 2, "test");
}

int fib(int n) {
  STACK();

  if (n > 1) {
    return fib(n - 1);
  } else {
    return test();
  }
}

auto main(int, char**) -> int {
  //

  STACK();

  fib(89);

  test();

  std::cout << "working\n";

  return 0;
}