// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <kineticpp/util.hpp>


namespace kineticpp::math {

template <typename T, typename JacStruct, typename JacLUStruct>
struct Eigen {

    template <size_t N>
    using Vector = ::Eigen::Matrix<T, N, 1>;

    using Scalar = T;
    using Jacobian = Vector<JacStruct::nnz>;

    class Solver {
        using SparseMatrix = ::Eigen::SparseMatrix<T>;
        using SparseLU = ::Eigen::SparseLU<SparseMatrix, ::Eigen::COLAMDOrdering<int>>;

        SparseMatrix A;
        SparseLU lu;

    public:
        Solver() : A(JacStruct::nrow, JacStruct::ncol) {
            std::array<::Eigen::Triplet<T>, JacStruct::nnz> triplets;
            JacStruct::for_ridx_cidx_rank([&](auto i, auto j, auto rank) {
                triplets[rank] = {i, j, 1.0};
            });
            A.setFromTriplets(triplets.begin(), triplets.end());
            lu.analyzePattern(A);
        }

        bool decompose(Jacobian &J) {
            std::array<::Eigen::Triplet<T>, JacStruct::nnz> triplets;
            JacStruct::for_ridx_cidx_rank([&](auto i, auto j, auto rank) {
                triplets[rank] = {i, j, J[rank]};
            });
            A.setFromTriplets(triplets.begin(), triplets.end());
            lu.factorize(A);
            return (lu.info() == ::Eigen::Success);
        }

        bool solve(Vector<JacLUStruct::nrow> &x) {
            x = lu.solve(x);
            return (lu.info() == ::Eigen::Success);
        }
    };

    template <auto N>
    static constexpr size_t size(Vector<N> &y) {
        return y.size();
    }

    // y <- 0
    template <auto N>
    static void zero(Vector<N> &y) {
        y.setZero();
    }

    // y <- x
    template <auto N>
    static void copy(Vector<N> &y, const Vector<N> &x, bool scrub = false) {
        if (scrub) {
            for (size_t i = 0; i < y.size(); ++i) {
                auto x_i = x(i);
                y(i) = util::is_nonzero(x_i) ? x_i : 0;
            }
        } else {
            y = x;
        }
    }

    // y <- alpha * y
    template <auto N>
    static void scale(Vector<N> &y, const Scalar alpha) {
        y *= alpha;
    }

    // y <- alpha * (y - x)
    template <auto N>
    static void aymx(Vector<N> &y, const Scalar alpha, const Vector<N> &x) {
        y = alpha * (y - x);
    }

    // y <- alpha * x + y
    template <auto N>
    static void axpy(Vector<N> &y, const Scalar alpha, const Vector<N> &x) {
        y = alpha * x + y;
    }

    // B <- I*alpha - A
    static void iama(Jacobian &B, const Scalar alpha, const Jacobian &A) {
        JacStruct::for_ridx_cidx_rank([&](auto i, auto j, auto rank) {
            if constexpr (i == j) {
                B[rank] = alpha - A[rank];
            } else {
                B[rank] = -A[rank];
            }
        });
    }
};

}  // namespace kineticpp::math
