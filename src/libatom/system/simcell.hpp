// #pragma once

// #include <cstddef>
// #include <cstdint>

// #include "bitsery/bitsery.h"
// #include "libatom/asserts.hpp"
// #include "libatom/system/atomvector.hpp"
// #include "libatom/system/simbox.hpp"
// #include "libatom/system/speciesmap.hpp"
// #include "libatom/utils.hpp"

// namespace otf {

//   /**
//    * @brief A SimCell is a collection of active/frozen atoms a species map and an OrthoSimBox
//    */
//   class SimCell {
//   public:
//     SpeciesMap map;

//     OrthoSimBox box;

//     AtomVector active;
//     AtomVector frozen;

//     /**
//      * @brief Construct a new Sim Cell object containing no atoms
//      */
//     explicit SimCell(OrthoSimBox const& box) : box{box} {}

//     /**
//      * @brief Get the total number of atoms in active + frozen
//      */
//     std::size_t size() const noexcept { return active.size() + frozen.size(); }

//     // VecN<double> active_disp(VecN<double> const &others) const;

//   protected:
//     friend class bitsery::Access;

//     /**
//      * @brief Construct a new SimCell object don't worry about class invariants, they will be
//      * restored in deserialization
//      */
//     SimCell() : box{OrthoSimBox::empty()} {}
//     /**
//      * @brief Bitsery serialisation
//      */
//     template <typename S> void serialize(S& s) {
//       ;
//       s(map);
//       s(box);
//       s(frozen);
//       s(active);
//     }
//   };

// }  // namespace otf