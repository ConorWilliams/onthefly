// #pragma once

// #include <fstream>

// #include "bitsery/adapter/stream.h"
// #include "bitsery/bitsery.h"
// #include "libatom/asserts.hpp"

// namespace otf {

//   /**
//    * @brief Dump a bitsery serialisable object to a file stream
//    *
//    * @code{.cpp}
//    *
//    * if (std::fstream f{"dump.bin", f.binary | f.trunc | f.out}; !f.is_open()) {
//    *    throw "Cannot open file for writing";
//    * } else {
//    *    otf::dump_binary(f, vec);
//    * }
//    *
//    * @endcode
//    */
//   template <typename Data> void dump_binary(std::fstream& file, Data const& data) {
//     ;

//     bitsery::Serializer<bitsery::OutputBufferedStreamAdapter> tmp{file};

//     tmp.object(data);

//     tmp.adapter().flush();
//   }

//   /**
//    * @brief Read a bitsery serialised object from a file
//    *
//    * @code{.cpp}
//    *
//    * if (std::fstream f{"dump.bin", s.binary | s.in}; !f.is_open()) {
//    *    throw "Cannot open file for reading";
//    * } else {
//    *    Obj obj;
//    *
//    *    while(true){
//    *        bool last = stream_binary(f, obj)
//    *
//    *        // Do something with each obj
//    *
//    *        if(last){
//    *            break
//    *        }
//    *    }
//    * }
//    *
//    * @endcode
//    *
//    * @return true If reached end of file
//    * @return false If did not read all of file (usefull for packed file)
//    */
//   template <typename Data> bool stream_binary(std::fstream& file, Data& data) {
//     ;

//     auto [err, done] = bitsery::quickDeserialization<bitsery::InputStreamAdapter>(file, data);

//     VERIFY(err == bitsery::ReaderError::NoError, "Deserialization failed");

//     return done;
//   }

// }  // namespace otf