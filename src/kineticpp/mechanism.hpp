// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <limits>
#include <type_traits>

#include "atom.hpp"
#include "expression.hpp"
#include "matrix.hpp"

namespace kineticpp {

template <typename T>
requires std::is_floating_point_v<T>
static constexpr bool is_nonzero(T elem) {
    return std::abs(elem) > (10 * std::numeric_limits<T>::epsilon());
}

template <typename T>
concept SpeciesId = std::is_enum_v<T>;

template <SpeciesId ID>
struct Species {
    using type = Species<ID>;
    using id_type = ID;

    const id_type id;
    const double mass;

    constexpr operator id_type() const { return static_cast<id_type>(id); }
};

template <typename L, typename R>
struct Equation {
    const L lhs;
    const R rhs;
};

template <typename L, typename R, typename Rate>
struct Reaction {
    const Equation<L, R> eqn;
    const Rate rate;
};

template <typename ID, Term T>
static constexpr auto operator||(ID lhs, T rhs) {
    return Species<ID> {lhs, atom::atomic_mass(rhs)};
}

// template <Term L, Term R>
template <typename L, typename R>
static constexpr auto operator>=(L lhs, R rhs) {
    return Equation<L, R> {lhs, rhs};
}

template <typename L, typename R, typename Rate>
static constexpr auto operator||(Equation<L, R> lhs, Rate rhs) {
    return Reaction<L, R, Rate> {lhs, rhs};
}

template <typename S, size_t N>
struct VariableSpecies : std::array<S, N> {};

template <typename T0, typename... TN>
VariableSpecies(T0, TN...) -> VariableSpecies<typename T0::type, 1 + sizeof...(TN)>;

template <typename S, size_t N>
struct FixedSpecies : std::array<S, N> {};

template <typename T0, typename... TN>
FixedSpecies(T0, TN...) -> FixedSpecies<typename T0::type, 1 + sizeof...(TN)>;


template <VariableSpecies Var, FixedSpecies Fix, auto... React>
class Mechanism {
public:
    using species_type = decltype(Var)::value_type;
    using species_id_type = species_type::id_type;

    static constexpr size_t nvar = Var.size();
    static constexpr size_t nfix = Fix.size();
    static constexpr size_t nspc = nvar + nfix;
    static constexpr size_t nrct = sizeof...(React);

    static constexpr bool is_var_spc(species_id_type id) {
        auto it = std::find(Var.begin(), Var.end(), id);
        return it != Var.end();
    }

    template <typename... Args>
    static constexpr auto rates(Args... args) {
        return std::array<double, nrct> {calc_rate(React.rate, args...)...};
    }

    static constexpr auto lhs_stoich() { return build_stoich_csr<true, false>(); }

    static constexpr auto rhs_stoich() { return build_stoich_csr<false, true>(); }

    static constexpr auto agg_stoich() { return build_stoich_csr<true, true>(); }

private:
    using term_type = std::pair<species_id_type, double>;

    template <Arithmetic Rate, typename... Args>
    static constexpr double calc_rate(const Rate &rconst, Args &&...) {
        return rconst;
    }

    template <typename Rate, typename... Args>
    requires std::is_invocable_v<Rate, Args...>
    static constexpr double calc_rate(const Rate &rfun, Args &&...args) {
        return rfun(args...);
    }

    template <typename Rate, typename... Args>
    static constexpr double calc_rate(const Rate &, Args &&...) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    static constexpr auto equation_terms(const species_id_type &expr) { return std::array {term_type {expr, 1.0}}; }

    template <Arithmetic A>
    static constexpr auto equation_terms(const Product<species_id_type, A> &expr) {
        return std::array {term_type {expr.lhs, expr.rhs}};
    }

    template <Arithmetic A>
    static constexpr auto equation_terms(const Product<A, species_id_type> &expr) {
        return std::array {term_type {expr.rhs, expr.lhs}};
    }

    template <Arithmetic A, typename RHS>
    static constexpr auto equation_terms(const Product<A, RHS> &expr) {
        double coef = expr.lhs;
        auto rhs_terms = equation_terms(expr.rhs);
        std::array<term_type, rhs_terms.size()> terms;
        for (size_t i = 0; i < terms.size(); ++i) {
            terms[i] = term_type {rhs_terms[i].first, coef * rhs_terms[i].second};
        }
        return terms;
    }

    template <typename LHS, Arithmetic A>
    static constexpr auto equation_terms(const Product<LHS, A> &expr) {
        auto lhs_terms = equation_terms(expr.lhs);
        double coef = expr.rhs;
        std::array<term_type, lhs_terms.size()> terms;
        for (size_t i = 0; i < terms.size(); ++i) {
            terms[i] = term_type {lhs_terms[i].first, coef * lhs_terms[i].second};
        }
        return terms;
    }

    template <typename LHS, typename RHS>
    static constexpr auto equation_terms(const Sum<LHS, RHS> &sum) {
        auto lhs_terms = equation_terms(sum.lhs);
        auto rhs_terms = equation_terms(sum.rhs);
        std::array<term_type, lhs_terms.size() + rhs_terms.size()> terms;
        std::copy(rhs_terms.begin(), rhs_terms.end(), std::copy(lhs_terms.begin(), lhs_terms.end(), terms.begin()));
        return terms;
    }

    template <bool lhs, bool rhs>
    static constexpr size_t count_nonzeros() {
        size_t nz = 0;
        foreach(
            [&](auto rct) {
                std::array<bool, nspc> row;
                row.fill(false);
                if constexpr (lhs && rhs) {
                    for (auto &term : equation_terms(rct.eqn.lhs)) {
                        if (is_var_spc(term.first)) {
                            row[spc_num(term.first)] = true;
                        }
                    }
                    for (auto &term : equation_terms(rct.eqn.rhs)) {
                        if (is_var_spc(term.first)) {
                            row[spc_num(term.first)] = true;
                        }
                    }
                } else if constexpr (lhs) {
                    for (auto &term : equation_terms(rct.eqn.lhs)) {
                        row[spc_num(term.first)] = true;
                    }
                } else if constexpr (rhs) {
                    for (auto &term : equation_terms(rct.eqn.rhs)) {
                        row[spc_num(term.first)] = true;
                    }
                }
                nz += std::count(row.begin(), row.end(), true);
            },
            React...);
        return nz;
    }

    template <bool lhs, bool rhs>
    static constexpr auto build_stoich_csr() {
        constexpr size_t rowsize = (lhs && rhs) ? nvar : nspc;
        constexpr size_t nz = count_nonzeros<lhs, rhs>();

        std::array<size_t, nrct + 1> ridx;
        std::array<size_t, nz> cols;
        std::array<double, nz> vals;

        size_t vals_idx = 0;
        size_t ridx_idx = 0;
        foreach(
            [&](auto rct) {
                std::array<double, rowsize> row;
                row.fill(0);
                if constexpr (lhs && rhs) {
                    for (auto &term : equation_terms(rct.eqn.lhs)) {
                        if (is_var_spc(term.first)) {
                            row[spc_num(term.first)] -= term.second;
                        }
                    }
                    for (auto &term : equation_terms(rct.eqn.rhs)) {
                        if (is_var_spc(term.first)) {
                            row[spc_num(term.first)] += term.second;
                        }
                    }
                } else if constexpr (lhs) {
                    for (auto &term : equation_terms(rct.eqn.lhs)) {
                        row[spc_num(term.first)] += term.second;
                    }
                } else if constexpr (rhs) {
                    for (auto &term : equation_terms(rct.eqn.rhs)) {
                        row[spc_num(term.first)] += term.second;
                    }
                }
                ridx[ridx_idx++] = vals_idx;
                for (size_t j = 0; j < row.size(); ++j) {
                    if (is_nonzero(row[j])) {
                        vals[vals_idx] = row[j];
                        cols[vals_idx] = j;
                        vals_idx++;
                    }
                }
            },
            React...);
        ridx[ridx_idx] = nz;
        return CsrMatrix<double, nrct, rowsize, nz> {ridx, cols, vals};
    }

    static constexpr size_t spc_num(species_id_type id) {
        auto it = std::find(spc_index.begin(), spc_index.end(), id);
        return it - spc_index.begin();
    }

    static constexpr auto declared_order() {
        std::array<species_id_type, nspc> order;
        for (size_t i = 0; i < Var.size(); ++i) {
            order[i] = static_cast<species_id_type>(Var[i]);
        }
        for (size_t i = 0; i < Fix.size(); ++i) {
            order[nvar + i] = static_cast<species_id_type>(Fix[i]);
        }
        return order;
    }

    static constexpr auto spc_index = declared_order();
};

}  // namespace kineticpp