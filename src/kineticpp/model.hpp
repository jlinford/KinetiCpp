// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "mechanism.hpp"
#include "solver/status.hpp"


namespace kineticpp {

template <typename T, typename Mech, template <typename, typename, typename, typename> typename S,
          template <typename, typename> typename LA>
class Model {
public:
    using MathLib = LA<T, typename Mech::jac_struct>;
    using VarConc = typename MathLib::Vector<Mech::nvar>;
    using FixConc = typename MathLib::Vector<Mech::nfix>;
    using Jacobian = typename MathLib::Jacobian;
    using Solver = S<MathLib, VarConc, FixConc, Jacobian>;

    static constexpr auto default_args = typename Solver::Parameters();

    static auto solve(VarConc &var, const FixConc &fix, double &h, const double t0, const double tend,
                      const Solver::Parameters &args = default_args) {
        return Solver::integrate(fun, jac, var, fix, h, t0, tend, args);
    }

    static auto solve(auto before_solve, VarConc &var, const FixConc &fix, double &h, const double t0,
                      const double tend, const double dt, const Solver::Parameters &args = default_args) {
        size_t step = 0;
        for (double t = t0; t < tend; t += dt, ++step) {
            before_solve(step, t, h, var, fix);
            auto errcode = Solver::integrate(fun, jac, var, fix, h, t, t + dt, args);
            if (errcode != solver::ErrorCode::success) {
                return errcode;
            }
        }
        return solver::ErrorCode::success;
    }

    static void fun(VarConc &du, const VarConc &var, const FixConc &fix, const double t) {
        using lhs = Mech::lhs_stoich;
        using agg = Mech::agg_stoich;
        auto rates = Mech::rates(t);

        std::array<double, Mech::nrct> rate_prod;
        lhs::for_ridx([&](auto i) {
            double prod = rates[i];
            lhs::for_cidx_val_in_row(i, [&](auto j, auto val) {
                if constexpr (j < Mech::nvar) {
                    prod *= std::pow(var[j], val);
                } else {
                    prod *= std::pow(fix[j - Mech::nvar], val);
                }
            });
            rate_prod[i] = prod;
        });

        MathLib::zero(du);
        agg::for_ridx_cidx_val([&](auto i, auto j, auto val) {
            du[j] += val * rate_prod[i];
        });
    }

    static void jac(Jacobian &J, const VarConc &var, const FixConc &fix, const double t) {
        using lhs = Mech::lhs_stoich;
        using agg = Mech::agg_stoich;
        using jac_struct = Mech::jac_struct;
        auto rates = Mech::rates(t);

        auto u = [&](auto j) {
            if constexpr (j < Mech::nvar) {
                return var[j];
            } else {
                return fix[j - Mech::nvar];
            }
        };

        std::array<double, lhs::nnz> B;
        size_t rank = 0;
        lhs::for_ridx_cidx([&](auto k, auto j) {
            double B_kj = rates[k];
            lhs::for_cidx_in_row(k, [&](auto jj) {
                constexpr double lhs_kjj = lhs::value(k, jj);
                if constexpr (jj == j) {
                    if constexpr (lhs_kjj == 2) {
                        B_kj *= 2 * u(jj);
                    } else if constexpr (lhs_kjj != 1) {
                        B_kj *= lhs_kjj * std::pow(u(jj), lhs_kjj - 1);
                    }
                } else {
                    if constexpr (lhs_kjj == 1) {
                        B_kj *= u(jj);
                    } else {
                        B_kj *= std::pow(u(jj), lhs_kjj);
                    }
                }
            });
            B[rank++] = B_kj;
        });

        rank = 0;
        jac_struct::for_ridx_cidx([&](auto i, auto j) {
            double sum = 0;
            for_constexpr<0, Mech::nrct>([&](auto k) {
                constexpr double A_ki = agg::value(k, i);
                if constexpr (is_nonzero(A_ki) && is_nonzero(lhs::value(k, j))) {
                    double B_kj = B[lhs::rank(k, j)];
                    sum += A_ki * B_kj;
                }
            });
            J[rank++] = sum;
        });
    }

};  // Model

}  // namespace kineticpp