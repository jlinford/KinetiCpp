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
        lhs::for_row_index([&](auto i) {
            double prod = rates[i];
            lhs::for_row_nz(i, [&](auto j, auto val) { prod *= std::pow(u[j], val); });
            rate_prod[i] = prod;
        });

        Solver::LinearAlgebra::zero(du);
        agg::for_nz([&](auto i, auto j, auto val) { du(j.value) += val * rate_prod[i]; });
    }

    static void jac(Jacobian &J, const Solution &u, const double t) {
        using lhs = Mech::lhs_stoich;
        using agg = Mech::agg_stoich;
        auto rates = Mech::rates(t);

        Solver::LinearAlgebra::zero(J);

        size_t rank = 0;
        lhs::for_nz([&](auto k, auto j, auto val) {
            double p = rates[k];
            lhs::for_row_nz(k, [&](auto ii, auto val) {
                if constexpr (ii == j) {
                    p *= val * std::pow(u[ii], val - 1);
                } else {
                    p *= std::pow(u[ii], val);
                }
            });
            for_constexpr<0, Mech::nvar>([&](auto i) {
                J(size_t(i),size_t(j)) += agg::value(k,i) * p;
            });
        });
        
    }

};  // Model

}  // namespace kineticpp