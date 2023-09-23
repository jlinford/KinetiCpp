// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>


namespace kineticpp {
namespace mathlib {

template <typename T, typename JS>
class EigenDense {
public:
    template <size_t N>
    using Vector = Eigen::Matrix<T, N, 1>;

    using Scalar = T;
    using JacStruct = JS;
    using Jacobian = Vector<JS::nnz>;

    struct Decomposition {
        using SparseMatrix = Eigen::SparseMatrix<T>;
        using SparseLU = Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>>;

        Decomposition() : A(JacStruct::nrow, JacStruct::ncol) {}

        void reset() { A.setZero(); }

        bool decompose(Jacobian &Anz) {
            std::array<Eigen::Triplet<T>, JacStruct::nnz> triplets;
            size_t rank = 0;
            JacStruct::for_ridx_cidx([&](auto i, auto j) {
                triplets[rank] = {i, j, Anz[rank]};
                ++rank;
            });
            A.setFromTriplets(triplets.begin(), triplets.end());
            solver.analyzePattern(A);
            solver.factorize(A);
            return (solver.info() == Eigen::Success);
        }

        auto solve(auto &x) {
            return solver.solve(x);
        }

        SparseMatrix A;
        SparseLU solver;
    };

public:
    static constexpr size_t size(auto &y) { return y.size(); }

    // y <- 0
    static void zero(auto &y) { y.setZero(); }

    // y <- x
    static void copy(auto &y, const auto &x) { y = x; }

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

    static bool decompose(Decomposition &decomp, Jacobian &Anz) { return decomp.decompose(Anz); }

    static void solve(Decomposition &decomp, auto &x) { x = decomp.solve(x); }
};


}  // namespace mathlib
}  // namespace kineticpp
