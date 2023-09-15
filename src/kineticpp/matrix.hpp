// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <type_traits>

#include "util.hpp"


namespace kineticpp {

template <typename T, size_t NR, size_t NC, size_t NNZ>
struct CsrMatrix {
    static constexpr size_t nrow = NR;
    static constexpr size_t ncol = NC;
    static constexpr size_t nnz = NNZ;

    const std::array<size_t, nrow + 1> ridx;
    const std::array<size_t, nnz> cols;
    const std::array<T, nnz> vals;

    constexpr T operator[](size_t row, size_t col) const {
        for (size_t ii = ridx[row]; ii < ridx[row + 1]; ++ii) {
            if (cols[ii] < col)
                continue;
            if (cols[ii] > col)
                return 0;
            return vals[ii];
        }
        return 0;
    }

    constexpr size_t rank(size_t row, size_t col) const {
        for (size_t ii = ridx[row]; ii < ridx[row + 1]; ++ii) {
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

    static constexpr size_t nrow = csr.nrow;
    static constexpr size_t ncol = csr.ncol;
    static constexpr size_t nnz = csr.nnz;

    static constexpr auto value(size_t row, size_t col) {
        return csr[row, col];
    }

    static constexpr size_t rank(size_t row, size_t col) {
        return csr.rank(row, col);
    }

    static constexpr void for_row_index(auto &&body) {
        for_constexpr<0, csr.nrow>([&](auto i) { body(std::integral_constant<size_t, i>()); });
    }

    static constexpr void for_row_col(auto &&body) {
        for_constexpr<0, csr.nrow>([&](auto i) {
            for_constexpr<csr.ridx[i], csr.ridx[i + 1]>([&](auto ii) {
                constexpr size_t j = csr.cols[ii];
                body(std::integral_constant<size_t, i>(), std::integral_constant<size_t, j>());
            });
        });
    }

    static constexpr void for_nz(auto &&body) {
        for_constexpr<0, csr.nrow>([&](auto i) {
            for_constexpr<csr.ridx[i], csr.ridx[i + 1]>([&](auto ii) {
                constexpr size_t j = csr.cols[ii];
                constexpr auto val = csr.vals[ii];
                body(std::integral_constant<size_t, i>(), std::integral_constant<size_t, j>(), val);
            });
        });
    }

    static constexpr void for_col(auto row, auto &&body) {
        for_constexpr<csr.ridx[row], csr.ridx[row + 1]>([&](auto ii) {
            constexpr size_t j = csr.cols[ii];
            body(std::integral_constant<size_t, j>());
        });
    }

    static constexpr void for_row_nz(auto row, auto &&body) {
        for_constexpr<csr.ridx[row], csr.ridx[row + 1]>([&](auto ii) {
            constexpr size_t j = csr.cols[ii];
            constexpr auto val = csr.vals[ii];
            body(std::integral_constant<size_t, j>(), val);
        });
    }
};


}  // namespace kineticpp
