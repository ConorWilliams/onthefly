#include "libatom/io/xyz.hpp"

#include <fmt/core.h>
#include <fmt/os.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <cstddef>
#include <optional>

#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::io {

  inline constexpr std::array symbols = {
      "00", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al",
      "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co",
      "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb",
      "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs",
      "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
      "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
      "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk",
      "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg",
  };

  void dump_xyz(fmt::ostream& file, SimCell const& atoms, std::string_view comment) {
    //

    // Write total number of atoms
    file.print("{}\n", atoms.size());

    // Build comment line

    VERIFY(std::find(comment.begin(), comment.end(), '\n') == comment.end(), "No newlines");

    file.print("{} ", comment);

    Vec3<floating> ext = atoms.extents();

    Vec3<bool> prd = atoms.periodic();

    file.print("Lattice=\"{} 0 0 0 {} 0 0 0 {}\" ", ext[0], ext[1], ext[2]);

    file.print("Periodic=\"{} {} {}\" ", prd[0], prd[1], prd[2]);

    file.print("Properties=atomic:I:1:species:S:1:pos:R:3:frozen:I:1\n");

    for (std::size_t i = 0; i < atoms.size(); i++) {
      //
      std::size_t n = atoms(AtomicNum{}, i);

      ASSERT(n < symbols.size(), "Atomic number out of bounds");

      file.print("{}\t{}\t{}\t{}\n", n, symbols[n], fmt::join(atoms(Position{}, i), "\t"),
                 (int)atoms(Frozen{}, i));
    }

    file.flush();
  }

  void stream_xyz(std::ifstream& file, SimCell& cell) {
    //
    std::size_t num_atoms;

    ASSERT(file, "Bad file supplied");

    file >> num_atoms;

    ASSERT(file && num_atoms == cell.size(), "Streaming the wrong number of atoms into the file");

    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // Skip to end of num line.
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // end of comments line.

    for (std::size_t i = 0; i < cell.size(); i++) {
      //
      file.ignore(1);  // Ignore species

      // Parse xyz
      file >> cell(Position{}, i)[0];
      file >> cell(Position{}, i)[1];
      file >> cell(Position{}, i)[2];

      ASSERT(file, "XYZ parsing error");

      file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // Skip to end of line.
    }
  }

}  // namespace otf::io