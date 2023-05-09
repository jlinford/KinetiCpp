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
        return Solver::integrate(fun, jac, u, h, t0, tend, args);
    }

    static auto solve(Solution &u, double &h, const double t0, const double tend, const double dt, auto before_solve,
                      const SolverParameters &args = default_args) {
        size_t step = 0;
        for (double t = t0; t < tend; t += dt, ++step) {
            before_solve(step, t, h, u);
            auto errcode = Model::solve(u, h, t, t + dt, args);
            if (errcode != solver::ErrorCode::success) {
                return errcode;
            }
        }
        return solver::ErrorCode::success;
    }

    static void fun(Solution &du, const Solution &u, const double t) {
        constexpr auto lhs = Mech::lhs_stoich();
        constexpr auto agg = Mech::agg_stoich();
        auto rates = Mech::rates(t);

        std::array<double, Mech::nrct> rate_prod;
        for_row<lhs>([&](auto i) {
            double prod = rates[i];
            for_nz<lhs, i>([&](auto j, auto val) { prod *= std::pow(u[j], val); });
            rate_prod[i] = prod;
        });

        Solver::LinearAlgebra::zero(du);
        for_nz<agg>([&](auto i, auto j, auto val) { du(j.value) += val * rate_prod[i]; });
    }

    static void jac(Jacobian &J, const Solution &u, const double t) {
        constexpr auto lhs = Mech::lhs_stoich();
        constexpr auto agg = Mech::agg_stoich();

        auto rates = Mech::rates(t);

        std::array<double, lhs.nnz> B;
        size_t rank = 0;
        for_nz<lhs>([&](auto i, auto j, auto val) {
            constexpr size_t jp1 = j + 1;
            double p = rates[i];
            // p *= prod(u[1:j-1].^lhs_stoich[1:j-1,i])
            for_constexpr<0, jp1 - 1>([&](auto k) {
                constexpr double lhs_exp = lhs[i, k];
                if constexpr (is_nonzero(lhs_exp)) {
                    p *= std::pow(u[k], lhs_exp);
                }
            });
            // p *= lhs_stoich[j,i] * u[j]^(lhs_stoich[j,i]-1)
            p *= val * std::pow(u[j], val - 1);
            // p *= prod(u[j+1:nspec].^lhs_stoich[j+1:nspec,i])
            for_constexpr<jp1, Mech::nspc>([&](auto k) {
                constexpr double lhs_exp = lhs[i, k];
                if constexpr (is_nonzero(lhs_exp)) {
                    p *= std::pow(u[k], lhs_exp);
                }
            });
            B[rank++] = p;
        });

        Solver::LinearAlgebra::zero(J);
        for_row_col<agg>([&](auto i, auto aj) {
            for_col<lhs, i>([&](auto lj) {
                // J[i,j] = sum(agg_stoich[:,i] .* B[:,j])
                if constexpr (lj < Mech::nvar) {
                    double sum = 0;
                    for_row<agg>([&](auto ii) {
                        constexpr size_t bi = lhs.rank(ii, lj);
                        if constexpr (bi < lhs.nnz) {
                            constexpr double agg_val = agg[ii, aj];
                            sum += agg_val * B[bi];
                        }
                    });
                    J(size_t(aj), size_t(lj)) = sum;
                }
            });
        });
    }

};  // Model

}  // namespace kineticpp