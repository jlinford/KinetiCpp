// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "mechanism.hpp"
#include "solver/status.hpp"


namespace kineticpp {

template <typename T, typename Mech, template <typename, typename, typename, typename> typename S,
          template <typename, typename> typename LA>
struct Model {

    using Scalar = T;
    using Timestep = double;
    using MathLib = LA<T, typename Mech::jac_struct>;
    using VarConc = typename MathLib::Vector<Mech::nvar>;
    using FixConc = typename MathLib::Vector<Mech::nfix>;
    using Jacobian = typename MathLib::Jacobian;
    using Solver = S<MathLib, VarConc, FixConc, Jacobian>;
    using SolverArgs = Solver::Args;

    static auto solve(VarConc &var, const FixConc &fix, const Timestep t0, const Timestep tend,
                      SolverArgs *args = NULL) {
        SolverArgs default_args;
        if (!args) {
            args = &default_args;
        } 
        return Solver::integrate(fun, jac, var, fix, t0, tend, *args);
    }

    static auto solve(auto before_solve, VarConc &var, const FixConc &fix, const Timestep t0, const Timestep tend,
                      const Timestep dt, SolverArgs *args = NULL) {
        SolverArgs default_args;
        size_t step = 0;
        if (!args) {
            args = &default_args;
        }
        for (Timestep t = t0; t < tend; t += dt, ++step) {
            before_solve(step, t, var, fix);
            auto errcode = solve(var, fix, t, t + dt, args);
            if (errcode != solver::ErrorCode::success) {
                return errcode;
            }
        }
        return solver::ErrorCode::success;
    }

    static void fun(VarConc &du, const VarConc &var, const FixConc &fix, const Timestep t) {
        using lhs = Mech::lhs_stoich;
        using agg = Mech::agg_stoich;
        auto rates = Mech::rates(t);

        std::array<Scalar, Mech::nrct> rate_prod;
        lhs::for_ridx([&](auto i) {
            Scalar prod = static_cast<Scalar>(rates[i]);
            lhs::for_cidx_val_in_row(i, [&](auto j, auto lhs_ij) {
                if constexpr (j < Mech::nvar) {
                    prod *= std::pow(var[j], static_cast<Scalar>(lhs_ij));
                } else {
                    prod *= std::pow(fix[j - Mech::nvar], static_cast<Scalar>(lhs_ij));
                }
            });
            rate_prod[i] = prod;
        });

        MathLib::zero(du);
        agg::for_ridx_cidx_val([&](auto i, auto j, auto agg_ij) {
            du[j] += agg_ij * rate_prod[i];
        });
    }

    static void jac(Jacobian &J, const VarConc &var, const FixConc &fix, const Timestep t) {
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

        std::array<Scalar, lhs::nnz> B;
        size_t rank = 0;
        lhs::for_ridx_cidx([&](auto k, auto j) {
            Scalar B_kj = static_cast<Scalar>(rates[k]);
            lhs::for_cidx_val_in_row(k, [&](auto jj, auto lhs_kjj) {
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
            Scalar sum = 0;
            for_constexpr<0, Mech::nrct>([&](auto k) {
                constexpr Scalar A_ki = agg::value(k, i);
                if constexpr (is_nonzero(A_ki) && is_nonzero(lhs::value(k, j))) {
                    Scalar B_kj = B[lhs::rank(k, j)];
                    sum += A_ki * B_kj;
                }
            });
            J[rank++] = sum;
        });
    }

};  // Model

}  // namespace kineticpp