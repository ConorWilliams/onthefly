#pragma once

#include <cstddef>
#include <fstream>
#include <string>

#include "local/classify.hpp"
#include "supercell.hpp"

class Stream {
  public:
    explicit Stream(Supercell const& init)
        : _cl(2.75),
          _t0(init),
          _raw{"olkmc.xyz", std::ios_base::app},
          _processed{"olkmc.txt", std::ios_base::app} {
        //
        _raw << std::setprecision(15);
        _processed << std::setprecision(15);
    }

    void operator()(Supercell const& cell,
                    int it,
                    double time,
                    double E0,
                    double bar,
                    double Ef,
                    double harm) {
        //

        double fe_sum = 0;
        double h_sum = 0;

        _cl(cell, _keys, _geos);

        _processed << it << ',' << time << ',' << E0 << ',' << bar << ',' << Ef << ',' << harm;

        for (std::size_t i = 0; i < cell.activ.size(); ++i) {
            switch (cell.activ[i].col.atomic) {
                case Colour::Fe:
                    fe_sum += norm_sq(cell.activ[i].vec - _t0.activ[i].vec);
                    break;
                case Colour::H:
                    h_sum += norm_sq(cell.activ[i].vec - _t0.activ[i].vec);

                    _processed << ',' << _keys[i].sdf[0];
                    _processed << ',' << cell.activ[i].vec[0];
                    _processed << ',' << cell.activ[i].vec[1];
                    _processed << ',' << cell.activ[i].vec[2];

                    break;
            }
        }

        //

        _processed << ',' << fe_sum << ',' << h_sum << '\n';
        _processed.flush();

        dump_raw(cell, it);
    }

    void dump_raw(Supercell const& cell, int it) {
        _raw << cell.size() << "\n";
        _raw << "It=" << it << " Lattice=\"" << cell.extents[0] << " 0.0 0.0 0.0 "
             << cell.extents[1] << " 0.0 0.0 0.0 " << cell.extents[2] << "\"";

        for (std::size_t i = 0; i < cell.activ.size(); ++i) {
            _raw << '\n' << cell.activ[i].col.atomic;
            _raw << ' ' << cell.activ[i].vec[0];
            _raw << ' ' << cell.activ[i].vec[1];
            _raw << ' ' << cell.activ[i].vec[2];
        }

        _raw << "\n";
        _raw.flush();
    }

  private:
    Classify _cl;

    std::vector<DiscreteKey> _keys;
    std::vector<Geometry> _geos;

    Supercell _t0;

    std::ofstream _raw;
    std::ofstream _processed;
};
