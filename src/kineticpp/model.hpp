// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "mechanism.hpp"
#include "solver/status.hpp"


namespace kineticpp {

template <typename Mech, template <size_t, typename> typename SolverT, typename LA>
class Model {
public:
    using Solver = SolverT<Mech::nspc, LA>;
    using SolverParameters = typename Solver::Parameters;
    using Solution = Solver::Vector;
    using Jacobian = Solver::Matrix;

    static constexpr auto default_args = SolverParameters();

    static auto solve(Solution &u, double &h, const double t0, const double tend,
                      const SolverParameters &args = default_args) {
        return Solver::integrate(fun, jac, u, h, t0, tend, default_args);
    }

    static auto solve(auto before_solve, Solution &u, double &h, const double t0, const double tend, const double dt,
                      const SolverParameters &args = default_args) {
        size_t step = 0;
        for (double t = t0; t < tend; t += dt, ++step) {
            before_solve(step, t, h, u);
            auto errcode = Solver::integrate(fun, jac, u, h, t, t + dt, args);
            if (errcode != solver::ErrorCode::success) {
                return errcode;
            }
        }
        return solver::ErrorCode::success;
    }

    static void fun(Solution &du, const Solution &u, const double t) {
        using lhs = Mech::lhs_stoich;
        using agg = Mech::agg_stoich;
        auto rates = Mech::rates(t);

        std::array<double, Mech::nrct> rate_prod;
        lhs::for_ridx([&](auto i) {
            double prod = rates[i];
            lhs::for_cidx_val_in_row(i, [&](auto j, auto val) { prod *= std::pow(u[j], val); });
            rate_prod[i] = prod;
        });

        Solver::LinearAlgebra::zero(du);
        agg::for_ridx_cidx_val([&](auto i, auto j, auto val) { du(j.value) += val * rate_prod[i]; });
    }

    static void jac(Jacobian &J, const Solution &u, const double t) {
        using lhs = Mech::lhs_stoich;
        using agg = Mech::agg_stoich;
        using jac_struct = Mech::jac_struct;
        auto rates = Mech::rates(t);

        Solver::LinearAlgebra::zero(J);

        std::array<double, Mech::nvar> B;
        size_t rank = 0;
        lhs::for_ridx_cidx([&](auto k, auto j) {
            double B_kj = rates[k];
            lhs::for_cidx_val_in_row(k, [&](auto jj, auto val) {
                if constexpr (jj == j) {
                    if constexpr (val == 2) {
                        B_kj *= 2 * u[jj];
                    } else if constexpr (val != 1) {
                        B_kj *= val * std::pow(u[jj], val - 1);
                    }
                } else {
                    if constexpr (val == 1) {
                        B_kj *= u[jj];
                    } else if constexpr (val == 2) {
                        B_kj *= u[jj]*u[jj];
                    } else {
                        B_kj *= std::pow(u[jj], val);
                    }
                }
            });
            B[rank++] = B_kj;
        });

        jac_struct::for_ridx_cidx([&](auto i, auto j) {
            double sum = 0;
            for_constexpr<0, Mech::nrct>([&](auto k) {
                constexpr double A_ki = agg::value(k, i);
                if constexpr (is_nonzero(A_ki) && is_nonzero(lhs::value(k, j))) {
                    double B_kj = B[lhs::rank(k, j)];
                    sum += A_ki * B_kj;
                }
            });
            J(size_t(i), size_t(j)) = sum;
        });
    }

};  // Model

}  // namespace kineticpp