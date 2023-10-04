// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <type_traits>

#include <kineticpp/util.hpp>


namespace kineticpp::math {

template <typename T, size_t NR, size_t NC, size_t NNZ>
struct CsrMatrix {
    using Index = size_t;
    using Value = T;

    static constexpr size_t nrow = NR;
    static constexpr size_t ncol = NC;
    static constexpr size_t ndiag = std::min(nrow, ncol);
    static constexpr size_t nnz = NNZ;

    const std::array<Index, nrow + 1> ridx;
    const std::array<Index, nnz> cols;
    const std::array<Index, ndiag> diag;
    const std::array<Value, nnz> vals;

    constexpr Value operator[](Index row, Index col) const {
        for (auto ii = ridx[row]; ii < ridx[row + 1]; ++ii) {
            if (cols[ii] < col)
                continue;
            if (cols[ii] > col)
                return 0;
            return vals[ii];
        }
        return 0;
    }

    constexpr Index rank(Index row, Index col) const {
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
    using Matrix = decltype(csr);
    using Index = Matrix::Index;
    using Value = Matrix::Value;

    template <auto x>
    using ConstexprIndex = std::integral_constant<Index, x>;

    template <auto x>
    using ConstexprValue = std::integral_constant<Value, x>;

    static constexpr size_t nrow = csr.nrow;
    static constexpr size_t ncol = csr.ncol;
    static constexpr size_t ndiag = csr.ndiag;
    static constexpr size_t nnz = csr.nnz;

    static constexpr Value value(auto row, auto col) { return csr[row, col]; }

    static constexpr Index rank(auto row, auto col) { return csr.rank(row, col); }

    static constexpr Index ridx(auto row) { return csr.ridx[row]; }

    static constexpr Index diag(auto row) { return csr.diag[row]; }

    static constexpr void for_ridx(auto &&body) {
        util::for_constexpr<0, csr.nrow>([&](auto i) {
            body(ConstexprIndex<i> {});
        });
    }

    static constexpr void for_ridx_reversed(auto &&body) {
        util::for_constexpr<csr.nrow, 0>([&](auto i) {
            body(ConstexprIndex<i> {});
        });
    }

    static constexpr void for_ridx_cidx(auto &&body) {
        util::for_constexpr<0, csr.nrow>([&](auto i) {
            util::for_constexpr<csr.ridx[i], csr.ridx[i + 1]>([&](auto ii) {
                body(ConstexprIndex<i> {}, ConstexprIndex<csr.cols[ii]> {});
            });
        });
    }

    static constexpr void for_ridx_cidx_val(auto &&body) {
        util::for_constexpr<0, csr.nrow>([&](auto i) {
            util::for_constexpr<csr.ridx[i], csr.ridx[i + 1]>([&](auto ii) {
                body(ConstexprIndex<i> {}, ConstexprIndex<csr.cols[ii]> {}, ConstexprValue<csr.vals[ii]> {});
            });
        });
    }

    static constexpr void for_ridx_cidx_rank(auto &&body) {
        util::for_constexpr<0, csr.nrow>([&](auto i) {
            util::for_constexpr<csr.ridx[i], csr.ridx[i + 1]>([&](auto ii) {
                body(ConstexprIndex<i> {}, ConstexprIndex<csr.cols[ii]> {}, ConstexprIndex<ii> {});
            });
        });
    }

    static constexpr void for_cidx_in_row(auto row, auto &&body) {
        util::for_constexpr<csr.ridx[row], csr.ridx[row + 1]>([&](auto ii) {
            body(ConstexprIndex<csr.cols[ii]> {});
        });
    }

    static constexpr void for_cidx_val_in_row(auto row, auto &&body) {
        util::for_constexpr<csr.ridx[row], csr.ridx[row + 1]>([&](auto ii) {
            body(ConstexprIndex<csr.cols[ii]> {}, ConstexprValue<csr.vals[ii]> {});
        });
    }

    static constexpr void for_cidx_rank_in_row(auto row, auto &&body) {
        util::for_constexpr<csr.ridx[row], csr.ridx[row + 1]>([&](auto ii) {
            body(ConstexprIndex<csr.cols[ii]> {}, ConstexprIndex<ii> {});
        });
    }

    static constexpr void for_cidx_rank_below_diag(auto row, auto &&body) {
        util::for_constexpr<csr.ridx[row], csr.diag[row]>([&](auto ii) {
            body(ConstexprIndex<csr.cols[ii]> {}, ConstexprIndex<ii> {});
        });
    }

    static constexpr void for_cidx_rank_above_diag(auto row, auto &&body) {
        util::for_constexpr<csr.diag[row] + 1, csr.ridx[row + 1]>([&](auto ii) {
            body(ConstexprIndex<csr.cols[ii]> {}, ConstexprIndex<ii> {});
        });
    }
};


}  // namespace kineticpp::math
