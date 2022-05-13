#include "utility.hpp"

#include <iostream>
#include <stdexcept>

// Error handler
[[noreturn]] __attribute__((noinline)) void assert_handler(int lineno,
                                                           char const *file,
                                                           char const *condition,
                                                           std::string message) {
    std::cerr << file << ':' << lineno << " assertion failed: " << condition << '\n';
    throw std::runtime_error(message);
}

// Throwing version of std::getline
std::istringstream safe_getline(std::ifstream &file) {
    if (std::string line; !std::getline(file, line)) {
        throw std::runtime_error("File terminated too soon");
    } else {
        return std::istringstream{line};
    }
}
