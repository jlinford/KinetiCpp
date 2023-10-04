// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <limits>
#include <type_traits>

#include <kineticpp/math/matrix.hpp>
#include <kineticpp/util.hpp>


namespace kineticpp::dsl {

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
    template <typename... Args>
    static constexpr auto rates(Args... args) {
        return std::array<double, nrct> {calc_rate(React.rate, args...)...};
    }

private:
    using species_type = decltype(Var)::value_type;
    using species_id_type = species_type::id_type;
    using term_type = std::pair<species_id_type, double>;

    static constexpr bool is_var_spc(species_id_type id) {
        auto it = std::find(Var.begin(), Var.end(), id);
        return it != Var.end();
    }

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
    static constexpr size_t count_stoich_nz() {
        size_t nz = 0;
        util::for_constexpr(
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
        constexpr size_t ncols = (lhs && rhs) ? nvar : nspc;
        constexpr size_t nz = count_stoich_nz<lhs, rhs>();
        constexpr size_t ndiag = std::min(nrct, ncols);

        std::array<size_t, nrct + 1> ridx;
        std::array<size_t, nz> cols;
        std::array<size_t, ndiag> diag;
        std::array<double, nz> vals;

        // Diagonal indices are unused in stoichiometric matrices
        std::fill(diag.begin(), diag.end(), nz);

        size_t vals_idx = 0;
        size_t ridx_idx = 0;
        util::for_constexpr(
            [&](auto rct) {
                std::array<double, ncols> row;
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
                ridx[ridx_idx] = vals_idx;
                for (size_t j = 0; j < row.size(); ++j) {
                    if (util::is_nonzero(row[j])) {
                        vals[vals_idx] = row[j];
                        cols[vals_idx] = j;
                        ++vals_idx;
                    }
                }
                ++ridx_idx;
            },
            React...);
        ridx[ridx_idx] = nz;
        return kineticpp::math::CsrMatrix<double, nrct, ncols, nz> {ridx, cols, diag, vals};
    }

    static constexpr size_t count_jac_nz() {
        size_t nz = 0;
        for (size_t i = 0; i < nvar; ++i) {
            for (size_t j = 0; j < nvar; ++j) {
                for (size_t k = 0; k < nrct; ++k) {
                    if (util::is_nonzero(agg_stoich::value(k, i)) && util::is_nonzero(lhs_stoich::value(k, j))) {
                        ++nz;
                        break;
                    }
                }
            }
        }
        return nz;
    }

    static constexpr auto build_jac_struct() {
        constexpr size_t nz = count_jac_nz();

        std::array<size_t, nvar + 1> ridx;
        std::array<size_t, nz> cols;
        std::array<size_t, nvar> diag;
        std::array<bool, nz> vals;

        size_t vals_idx = 0;
        size_t ridx_idx = 0;
        for (size_t i = 0; i < nvar; ++i) {
            std::array<bool, nvar> row;
            row.fill(false);
            for (size_t j = 0; j < nvar; ++j) {
                if (i == j) {
                    row[j] = true;
                } else {
                    for (size_t k = 0; k < nrct; ++k) {
                        if (util::is_nonzero(agg_stoich::value(k, i)) && util::is_nonzero(lhs_stoich::value(k, j))) {
                            row[j] = true;
                            break;
                        }
                    }
                }
            }
            ridx[ridx_idx] = vals_idx;
            for (size_t j = 0; j < row.size(); ++j) {
                if (i == j) {
                    vals[vals_idx] = row[j];
                    cols[vals_idx] = j;
                    diag[ridx_idx] = vals_idx;
                    ++vals_idx;
                } else if (row[j]) {
                    vals[vals_idx] = row[j];
                    cols[vals_idx] = j;
                    ++vals_idx;
                }
            }
            ++ridx_idx;
        }
        ridx[ridx_idx] = nz;
        return kineticpp::math::CsrMatrix<bool, nvar, nvar, nz> {ridx, cols, diag, vals};
    }

    static constexpr auto dense_jac_lu_struct() {
        std::array<bool, nvar * nvar> dense;
        dense.fill(false);

        jac_struct::for_ridx_cidx([&](auto i, auto j) {
            dense[i * nvar + j] = true;
        });

        for (size_t j = 0; j < nvar - 1; ++j) {
            for (size_t i = j + 1; i < nvar; ++i) {
                if (dense[i * nvar + j]) {
                    for (size_t k = j; k < nvar; ++k) {
                        if (dense[j * nvar + k]) {
                            dense[i * nvar + k] = true;
                        }
                    }
                }
            }
        }
        return dense;
    }

    static constexpr auto build_jac_lu_struct() {
        constexpr auto dense = dense_jac_lu_struct();
        constexpr size_t nz = std::count(dense.begin(), dense.end(), true);

        std::array<size_t, nvar + 1> ridx;
        std::array<size_t, nz> cols;
        std::array<size_t, nvar> diag;
        std::array<bool, nz> vals;

        size_t vals_idx = 0;
        size_t ridx_idx = 0;
        for (size_t i = 0; i < nvar; ++i) {
            ridx[ridx_idx] = vals_idx;
            for (size_t j = 0; j < nvar; ++j) {
                if (i == j) {
                    vals[vals_idx] = dense[i * nvar + j];
                    cols[vals_idx] = j;
                    diag[ridx_idx] = vals_idx;
                    ++vals_idx;
                } else if (dense[i * nvar + j]) {
                    vals[vals_idx] = dense[i * nvar + j];
                    cols[vals_idx] = j;
                    ++vals_idx;
                }
            }
            ++ridx_idx;
        }
        ridx[ridx_idx] = nz;
        return kineticpp::math::CsrMatrix<bool, nvar, nvar, nz> {ridx, cols, diag, vals};
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

public:
    using lhs_stoich = kineticpp::math::ConstexprCsrMatrix<build_stoich_csr<true, false>()>;
    using rhs_stoich = kineticpp::math::ConstexprCsrMatrix<build_stoich_csr<false, true>()>;
    using agg_stoich = kineticpp::math::ConstexprCsrMatrix<build_stoich_csr<true, true>()>;
    using jac_struct = kineticpp::math::ConstexprCsrMatrix<build_jac_struct()>;
    using jac_lu_struct = kineticpp::math::ConstexprCsrMatrix<build_jac_lu_struct()>;

    static constexpr size_t nvar = Var.size();
    static constexpr size_t nfix = Fix.size();
    static constexpr size_t nspc = nvar + nfix;
    static constexpr size_t nrct = sizeof...(React);
    static constexpr size_t njac = jac_struct::nnz;
    static constexpr size_t njac_lu = jac_lu_struct::nnz;
};

}  // namespace kineticpp::dsl
