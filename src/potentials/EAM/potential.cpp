

#include "potentials/EAM/potential.hpp"

#include <math.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include "config.hpp"
#include "potentials/EAM/data.hpp"
#include "potentials/neigh_reduce.hpp"
#include "supercell.hpp"

PotentialEAM::PotentialEAM(TabEAM &&data) : _data(std::make_shared<TabEAM>(std::move(data))) {}

std::unique_ptr<PotentialBase> PotentialEAM::clone() const {
    return std::make_unique<PotentialEAM>(*this);
}

// Compute energy
double PotentialEAM::energy(Supercell const &cell) {
    _reduce_energy.load(rcut(), cell);

    double v_sum = 0;
    double f_sum = 0;

    for (auto a = _reduce_energy.begin_activ(); a != _reduce_energy.begin_bound(); ++a) {
        double rho = 0;

        _reduce_energy.neigh_reduce(a, [&](auto b, double r, Vec3<double> const &) {
            v_sum += _data->v(a->col.atomic, b->col.atomic)(r);
            rho += _data->phi(b->col.atomic, a->col.atomic)(r);
        });

        f_sum += _data->f(a->col.atomic)(rho);
    }

    return (0.5 * v_sum) + f_sum;
}

// Compute force
void PotentialEAM::gradient(Supercell const &cell, VecN<double> &out) {
    _reduce_grad.load(rcut(), cell);

    // First sum computes density  at each atom, runs over active+boundary atoms
    for (auto b = _reduce_grad.begin_activ(); b != _reduce_grad.begin_ghost(); ++b) {
        double rho = 0;

        // Computes rho at atom
        _reduce_grad.neigh_reduce(
            b, [&](auto a, double r) { rho += _data->phi(a->col.atomic, b->col.atomic)(r); });

        // Compute F'(rho) at atom
        b->fp_rho = _data->f(b->col.atomic).grad(rho);
    }

    // Update density of ghost atoms
    _reduce_grad.broadcast_ghost_data();

    out.resize(3 * cell.activ.size());  // Usually a no-op

    // Counter must be declared outside loop
    std::size_t i = 0;
    // Second sum computes force, only runs over active atoms
    for (auto g = _reduce_grad.begin_activ(); g != _reduce_grad.begin_bound(); ++g) {
        out[3 * i + 0] = 0;
        out[3 * i + 1] = 0;
        out[3 * i + 2] = 0;

        // Finds R^{\alpha\gamma}
        _reduce_grad.neigh_reduce(g, [&](auto a, double r, Vec3<double> const &dr) {
            double mag = _data->v(a->col.atomic, g->col.atomic).grad(r) +

                         g->fp_rho * _data->phi(a->col.atomic, g->col.atomic).grad(r) +

                         a->fp_rho * _data->phi(g->col.atomic, a->col.atomic).grad(r);

            mag /= r;

            out[3 * i + 0] += mag * dr[0];
            out[3 * i + 1] += mag * dr[1];
            out[3 * i + 2] += mag * dr[2];
        });

        ++i;
    }
}

void PotentialEAM::hessian(Supercell const &cell, MatN<double> &out) {
    _reduce_hess.load(rcut(), cell);

    std::size_t bc = 0;
    // First sum computes density  at each atom, runs over active+boundary atoms
    for (auto b = _reduce_hess.begin_activ(); b != _reduce_hess.begin_ghost(); ++b, ++bc) {
        // Computes rho and mu at atom and tags index
        _reduce_hess.neigh_reduce(b, [&](auto a, double r, Vec3<double> const &dr) {
            //
            b->rho += _data->phi(a->col.atomic, b->col.atomic)(r);
            b->mu += _data->phi(a->col.atomic, b->col.atomic).grad(r) * dr / r;  // A.13
            b->idx = bc;
        });
    }

    // Update rho & mu of ghost atoms, additionally each ghost contains index of idx
    _reduce_hess.broadcast_ghost_data();

    out = MatN<double>::Zero(3 * cell.activ.size(), 3 * cell.activ.size());

    // Require double neighbour-loop for overlap terms; accelerate with cache
    std::vector<NeighReduce<Hess>::neigh_atom *> neigh_cache;

    // Second sums computes hessian, only runs over active atoms
    for (auto x = _reduce_hess.begin_activ(); x != _reduce_hess.begin_bound(); ++x) {
        // Block-diagonal hessian-elements x ~ gamma
        double ddFg = _data->f(x->col.atomic).grad2(x->rho);
        // A.15 pre sum term
        for (size_t j = 0; j < 3; j++) {
            for (size_t i = 0; i < 3; i++) {
                out(3 * x->idx + i, 3 * x->idx + j) = ddFg * x->mu(i) * x->mu(j);
            }
        }

        neigh_cache.clear();

        // Note: dr = r^{ax}
        _reduce_hess.neigh_reduce(x, [&](auto a, double r, Vec3<double> const &dr) {
            // First compute sum over neigh for block-diagonal hessian-elements

            // A.14
            double A = _data->v(a->col.atomic, x->col.atomic).grad(r) +

                       _data->f(x->col.atomic).grad(x->rho)
                           * _data->phi(a->col.atomic, x->col.atomic).grad(r)
                       +

                       _data->f(a->col.atomic).grad(a->rho)
                           * _data->phi(x->col.atomic, a->col.atomic).grad(r);

            A /= r;

            // A.14
            double B = _data->v(a->col.atomic, x->col.atomic).grad2(r) +

                       _data->f(x->col.atomic).grad(x->rho)
                           * _data->phi(a->col.atomic, x->col.atomic).grad2(r)
                       +

                       _data->f(a->col.atomic).grad(a->rho)
                           * _data->phi(x->col.atomic, a->col.atomic).grad2(r);

            double s = A - B
                       - _data->f(a->col.atomic).grad2(a->rho)
                             * _data->phi(x->col.atomic, a->col.atomic).grad(r)
                             * _data->phi(x->col.atomic, a->col.atomic).grad(r);

            s /= (r * r);

            // A.15 inside sum
            for (size_t j = 0; j < 3; j++) {
                // Isotropic component
                out(3 * x->idx + j, 3 * x->idx + j) += A;

                for (size_t i = 0; i < 3; i++) {
                    // Sign flip for r^{xa}
                    out(3 * x->idx + i, 3 * x->idx + j) -= s * dr(i) * dr(j);
                }
            }

            /////////////////////////////////////////////////////////////

            if (a->idx < cell.activ.size()) {
                // Cache neighbours inside active
                neigh_cache.push_back(a);

                // Now compute \eta \gamma terms for r^{\eta\gamma} < r_cut

                // let:   \eta \get x
                // let: \gamma \get a

                // Therefore : dr = r^{\gamma\eta}

                double s2 = (B - A) / (r * r);

                double ur = _data->f(x->col.atomic).grad2(x->rho)
                            * _data->phi(a->col.atomic, x->col.atomic).grad(r) / r;

                double ru = _data->f(a->col.atomic).grad2(a->rho)
                            * _data->phi(x->col.atomic, a->col.atomic).grad(r) / r;

                for (size_t j = 0; j < 3; j++) {
                    // Isotropic component
                    out(3 * x->idx + j, 3 * a->idx + j) -= A;

                    for (size_t i = 0; i < 3; i++) {
                        out(3 * x->idx + i, 3 * a->idx + j) -= s2 * dr(i) * dr(j);

                        out(3 * x->idx + i, 3 * a->idx + j) -= ur * x->mu(i) * dr(j);
                        out(3 * x->idx + i, 3 * a->idx + j) += ru * dr(i) * a->mu(j);
                    }
                }
            }
        });

        // Compute overlap terms between pairs coupled by x ~ \alpha inside activ

        for (auto eta : neigh_cache) {
            for (auto gam : neigh_cache) {
                if (eta != gam) {
                    Vec3<double> rna = x->vec - eta->vec;
                    Vec3<double> rag = gam->vec - x->vec;

                    double rn = norm(rna);
                    double rg = norm(rag);

                    double mag = ddFg / (rn * rg)
                                 * _data->phi(gam->col.atomic, x->col.atomic).grad(rg)
                                 * _data->phi(eta->col.atomic, x->col.atomic).grad(rn);

                    for (size_t j = 0; j < 3; j++) {
                        for (size_t i = 0; i < 3; i++) {
                            out(3 * eta->idx + i, 3 * gam->idx + j) -= mag * rna(i) * rag(j);
                        }
                    }
                }
            }
        }
    }

    // Third sum computes overlap term between pairs coupled by \alpha outside activ (inside bound)
    for (auto a = _reduce_hess.begin_bound(); a != _reduce_hess.begin_ghost(); ++a) {
        neigh_cache.clear();

        _reduce_hess.neigh_reduce(a, [&](auto n) {
            if (n->idx < cell.activ.size()) {
                neigh_cache.push_back(n);
            }
        });

        double ddFa = _data->f(a->col.atomic).grad2(a->rho);

        for (auto eta : neigh_cache) {
            for (auto gam : neigh_cache) {
                if (eta != gam) {
                    Vec3<double> rna = a->vec - eta->vec;
                    Vec3<double> rag = gam->vec - a->vec;

                    double rn = norm(rna);
                    double rg = norm(rag);

                    double mag = ddFa / (rn * rg)
                                 * _data->phi(gam->col.atomic, a->col.atomic).grad(rg)
                                 * _data->phi(eta->col.atomic, a->col.atomic).grad(rn);

                    for (size_t j = 0; j < 3; j++) {
                        for (size_t i = 0; i < 3; i++) {
                            out(3 * eta->idx + i, 3 * gam->idx + j) -= mag * rna(i) * rag(j);
                        }
                    }
                }
            }
        }
    }

    // Transform into mass-weighted hessian
    for (auto i = _reduce_hess.begin_activ(); i != _reduce_hess.begin_bound(); ++i) {
        for (auto j = _reduce_hess.begin_activ(); j != _reduce_hess.begin_bound(); ++j) {
            out.block<3, 3>(3 * j->idx, 3 * i->idx)
                *= 1 / std::sqrt(_data->mass(i->col.atomic) * _data->mass(j->col.atomic));
        }
    }
}
