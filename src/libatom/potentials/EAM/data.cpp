#include "libatom/potentials/EAM/data.hpp"

#include <cstddef>
#include <fstream>
#include <limits>

#include "libatom/asserts.hpp"
#include "libatom/utils.hpp"

namespace otf::potentials {

  /**
   * @brief Throwing version of std::getline
   */
  static std::istringstream getline(std::ifstream& in) {
    if (std::string line; !std::getline(in, line)) {
      throw std::runtime_error("File terminated too soon");
    } else {
      return std::istringstream{line};
    }
  }

  /**
   * @brief Read N elements in lines of length K into a vector
   */
  static std::vector<floating> read_chunked(std::ifstream& file, std::size_t n, std::size_t k) {
    //
    std::vector<floating> raw;

    std::istringstream line;

    for (std::size_t i = 0; i < n; ++i) {
      if (i % k == 0) {
        line = getline(file);
      }

      floating tmp;
      line >> tmp;
      raw.push_back(tmp);
    }

    return raw;
  }

  DataEAM::DataEAM(std::ifstream in) {
    //
    VERIFY(in.good(), "Could not open eam in");

    m_atomic2idx.fill(std::numeric_limits<std::size_t>::max());
    m_atomic2mass.fill(0);

    // Skip to 4th line
    for (int i = 0; i < 3; ++i) {
      getline(in);
    }

    // Parse number of species
    getline(in) >> m_num_species;

    // Allocate space for splines
    m_f.resize(m_num_species);
    m_phi.resize(m_num_species * m_num_species);
    m_v.resize(m_num_species * m_num_species);

    // Temporaries
    std::size_t numP, numR;
    floating delP, delR;

    // Parse tabulation info
    getline(in) >> numP >> delP >> numR >> delR >> m_rcut;

    for (std::size_t i = 0; i < m_num_species; ++i) {
      {  // Read species info
        std::size_t atomic;
        floating mass;
        getline(in) >> atomic >> mass;

        VERIFY(atomic < 112, "Not a valid atomic number");

        m_atomic2idx[atomic] = i;
        m_atomic2mass[atomic] = mass;
      }

      // Read F
      m_f[i] = Spline{read_chunked(in, numP, 5), delP};

      // Read phi
      for (std::size_t j = 0; j < m_num_species; ++j) {
        m_phi[index(i, j)] = Spline{read_chunked(in, numR, 5), delR};
      }
    }

    // Read v ***IMPORTANT**** tabulated as r*v and symmetric

    for (std::size_t i = 0; i < m_num_species; ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        //
        std::vector raw = read_chunked(in, numP, 5);

        for (std::size_t k = 0; k < raw.size(); k++) {
          raw[k] /= delR * k;
        }

        VERIFY(raw.size() > 0, "no elements!");

        raw[0] = raw[1];  // Fixup divide by zero for spline

        m_v[sym_index(i, j)] = Spline{std::move(raw), delR};
      }
    }
  }

}  // namespace otf::potentials