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

  static void dump(fmt::ostream& file, SpeciesMap const& map, AtomVector const& atoms) {
    STACK();

    for (size_t i = 0; i < atoms.size(); i++) {
      file.print("{}\t{}\n", map.z2species(atoms.z()[i]), fmt::join(atoms.x().col(i), "\t"));
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
      file.print("Lattice=\"{} 0 0 0 {} 0 0 0 {}\" ", cell.box.extents()[0], cell.box.extents()[1],
                 cell.box.extents()[2]);
    } else if constexpr (spatial_dims == 2) {
      file.print("Lattice=\"{} 0 0 {}\" ", cell.box.extents()[0], cell.box.extents()[1]);
    }

    file.print("Properties=species:S:1:pos:R:{}\n", spatial_dims);

    dump(file, cell.map, cell.active);
    dump(file, cell.map, cell.frozen);

    file.flush();
  }

}  // namespace otf