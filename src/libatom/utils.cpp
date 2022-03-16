
#include <string_view>

namespace otf {

  std::string_view common_prefix(std::string_view a, std::string_view b) {
    for (std::size_t i = 0; i < std::min(a.size(), b.size()); i++) {
      if (a[i] != b[i]) {
        return a.substr(0, i);
      }
    }
    return a.substr(0, std::min(a.size(), b.size()));
  }

}  // namespace otf
