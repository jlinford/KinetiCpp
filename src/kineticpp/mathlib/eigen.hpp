// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <kineticpp/util.hpp>

namespace kineticpp::mathlib {

template <typename T, typename JS>
struct Eigen {

    template <size_t N>
    using Vector = ::Eigen::Matrix<T, N, 1>;

    using Scalar = T;
    using JacStruct = JS;
    using Jacobian = Vector<JS::nnz>;

    class Solver {
        using SparseMatrix = ::Eigen::SparseMatrix<T>;
        using SparseLU = ::Eigen::SparseLU<SparseMatrix, ::Eigen::COLAMDOrdering<int>>;

        SparseMatrix A;
        SparseLU lu;

    public:
        Solver() : A(JacStruct::nrow, JacStruct::ncol) {}

        auto decompose(Jacobian &Anz) {
            std::array<::Eigen::Triplet<T>, JacStruct::nnz> triplets;
            size_t rank = 0;
            JacStruct::for_ridx_cidx([&](auto i, auto j) {
                triplets[rank] = {i, j, Anz[rank]};
                ++rank;
            });
            A.setFromTriplets(triplets.begin(), triplets.end());
            lu.analyzePattern(A);
            lu.factorize(A);
            return (lu.info() == ::Eigen::Success);
        }

        auto solve(auto &x) {
            x = lu.solve(x);
            return (lu.info() == ::Eigen::Success);
        }
    };

    static constexpr size_t size(auto &y) { return y.size(); }

    // y <- 0
    static void zero(auto &y) { y.setZero(); }

    // y <- x
    static void copy(auto &y, const auto &x, bool scrub = false) {
        if (scrub) {
            for (size_t i = 0; i < y.size(); ++i) {
                auto x_i = x(i);
                y(i) = is_nonzero(x_i) ? x_i : 0;
            }
        } else {
            y = x;
        }
    }

    // y <- alpha * y
    static void scale(auto &y, const auto alpha) { y *= alpha; }

    // y <- alpha * (y - x)
    static void aymx(auto &y, const auto alpha, const auto &x) { y = alpha * (y - x); }

    // y <- alpha * x + y
    static void axpy(auto &y, const auto alpha, const auto &x) { y = alpha * x + y; }

    // B <- I*alpha - A
    static void iama(Jacobian &B, auto const alpha, const Jacobian &A) {
        size_t rank = 0;
        JacStruct::for_ridx_cidx([&](auto i, auto j) {
            if constexpr (i == j) {
                B[rank] = alpha - A[rank];
            } else {
                B[rank] = -A[rank];
            }
            ++rank;
        });
    }
};

}  // namespace kineticpp::mathlib
