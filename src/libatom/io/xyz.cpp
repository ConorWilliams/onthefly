#include "libatom/io/xyz.hpp"

#include <fmt/core.h>

#include <optional>

#include "fmt/os.h"
#include "fmt/ranges.h"
#include "libatom/asserts.hpp"
#include "libatom/system/atomvector.hpp"
#include "libatom/system/simcell.hpp"
#include "libatom/utils.hpp"

namespace otf {

  static void dump(fmt::ostream& file, AtomVector const& atoms) {
    STACK();

    for (size_t i = 0; i < atoms.size(); i++) {
      if (Symbol s = atoms.z2species(atoms.z()[i]); s[1]) {
        file.print("{}{}\t{}\n", s[0], s[1], fmt::join(atoms.x().col(i), "\t"));
      } else {
        file.print("{}\t{}\n", s[0], fmt::join(atoms.x().col(i), "\t"));
      }
    }
  }

  void dump_xyz(fmt::ostream& file, SimCell const& cell, std::optional<flt_t> time) {
    //

    STACK();

    // Write total number of atoms
    file.print("{}\n", cell.size());

    // Build comment line

    if (time) {
      file.print("Time={} ", *time);
    }

    if constexpr (spatial_dims == 3) {
      file.print("Lattice=\"{} 0 0 0 {} 0 0 0 {}\" ", cell.extents()[0], cell.extents()[1],
                 cell.extents()[2]);
    } else if constexpr (spatial_dims == 2) {
      file.print("Lattice=\"{} 0 0 {}\" ", cell.extents()[0], cell.extents()[1]);
    }

    file.print("Properties=species:S:1:pos:R:{}\n", spatial_dims);

    dump(file, cell.active());
    dump(file, cell.frozen());

    file.flush();
  }

}  // namespace otf