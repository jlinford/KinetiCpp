// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <type_traits>

#include "util.hpp"


namespace kineticpp {

template <typename T, size_t NR, size_t NC, size_t NNZ>
struct CsrMatrix {
    using index_type = size_t;
    using value_type = T;

    static constexpr size_t nrow = NR;
    static constexpr size_t ncol = NC;
    static constexpr size_t ndiag = std::min(nrow, ncol);
    static constexpr size_t nnz = NNZ;

    const std::array<index_type, nrow + 1> ridx;
    const std::array<index_type, nnz> cols;
    const std::array<index_type, ndiag> diag;
    const std::array<value_type, nnz> vals;

    constexpr value_type operator[](index_type row, index_type col) const {
        for (auto ii = ridx[row]; ii < ridx[row + 1]; ++ii) {
            if (cols[ii] < col)
                continue;
            if (cols[ii] > col)
                return 0;
            return vals[ii];
        }
        return 0;
    }

    constexpr index_type rank(index_type row, index_type col) const {
        for (auto ii = ridx[row]; ii < ridx[row + 1]; ++ii) {
            if (cols[ii] < col)
                continue;
            if (cols[ii] > col)
                return nnz;
            return ii;
        }
        return nnz;
    }
};

template <CsrMatrix csr>
struct ConstexprCsrMatrix {
    using csr_type = decltype(csr);
    using index_type = csr_type::index_type;
    using value_type = csr_type::value_type;

    template <auto x>
    using constexpr_index = std::integral_constant<index_type, x>;

    template <auto x>
    using constexpr_value = std::integral_constant<value_type, x>;

    static constexpr size_t nrow = csr.nrow;
    static constexpr size_t ncol = csr.ncol;
    static constexpr size_t ndiag = csr.ndiag;
    static constexpr size_t nnz = csr.nnz;

    static constexpr value_type value(auto row, auto col) { return csr[row, col]; }

    static constexpr index_type rank(auto row, auto col) { return csr.rank(row, col); }

    static constexpr index_type ridx(auto row) { return csr.ridx[row]; }

    static constexpr index_type diag(auto row) { return csr.diag[row]; }

    static constexpr void for_ridx(auto &&body) {
        for_constexpr<0, csr.nrow>([&](auto i) {
            body(constexpr_index<i> {});
        });
    }

    static constexpr void for_ridx_reversed(auto &&body) {
        for_constexpr<csr.nrow, 0>([&](auto i) {
            body(constexpr_index<i> {});
        });
    }

    static constexpr void for_ridx_cidx(auto &&body) {
        for_constexpr<0, csr.nrow>([&](auto i) {
            for_constexpr<csr.ridx[i], csr.ridx[i + 1]>([&](auto ii) {
                body(constexpr_index<i> {}, constexpr_index<csr.cols[ii]> {});
            });
        });
    }

    static constexpr void for_ridx_cidx_val(auto &&body) {
        for_constexpr<0, csr.nrow>([&](auto i) {
            for_constexpr<csr.ridx[i], csr.ridx[i + 1]>([&](auto ii) {
                body(constexpr_index<i> {}, constexpr_index<csr.cols[ii]> {}, constexpr_value<csr.vals[ii]> {});
            });
        });
    }

    static constexpr void for_ridx_cidx_rank(auto &&body) {
        for_constexpr<0, csr.nrow>([&](auto i) {
            for_constexpr<csr.ridx[i], csr.ridx[i + 1]>([&](auto ii) {
                body(constexpr_index<i> {}, constexpr_index<csr.cols[ii]> {}, constexpr_index<ii> {});
            });
        });
    }

    static constexpr void for_cidx_in_row(auto row, auto &&body) {
        for_constexpr<csr.ridx[row], csr.ridx[row + 1]>([&](auto ii) {
            body(constexpr_index<csr.cols[ii]> {});
        });
    }

    static constexpr void for_cidx_val_in_row(auto row, auto &&body) {
        for_constexpr<csr.ridx[row], csr.ridx[row + 1]>([&](auto ii) {
            body(constexpr_index<csr.cols[ii]> {}, constexpr_value<csr.vals[ii]> {});
        });
    }

    static constexpr void for_cidx_rank_in_row(auto row, auto &&body) {
        for_constexpr<csr.ridx[row], csr.ridx[row + 1]>([&](auto ii) {
            body(constexpr_index<csr.cols[ii]> {}, constexpr_index<ii> {});
        });
    }

    static constexpr void for_cidx_rank_below_diag(auto row, auto &&body) {
        for_constexpr<csr.ridx[row], csr.diag[row]>([&](auto ii) {
            body(constexpr_index<csr.cols[ii]> {}, constexpr_index<ii> {});
        });
    }

    static constexpr void for_cidx_rank_above_diag(auto row, auto &&body) {
        for_constexpr<csr.diag[row] + 1, csr.ridx[row + 1]>([&](auto ii) {
            body(constexpr_index<csr.cols[ii]> {}, constexpr_index<ii> {});
        });
    }
};


}  // namespace kineticpp
