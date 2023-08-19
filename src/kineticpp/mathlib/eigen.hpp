// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <Eigen/Dense>


namespace kineticpp {
namespace mathlib {


template <typename T>
struct EigenDense {
    using Scalar = T;

    template <size_t N>
    using Vector = Eigen::Matrix<T, N, 1>;

    template <size_t Rows, size_t Cols>
    using Matrix = Eigen::Matrix<T, Rows, Cols>;

    template <size_t Rows>
    static constexpr size_t size(Vector<Rows> &y) {
        return y.size();
    }

    // y <- 0
    static void zero(auto &y) { y.setZero(); }

    // y <- x
    static void copy(auto &y, const auto &x) { y = x; }

    // y <- alpha * y
    static void scale(auto &y, const Scalar alpha) { y *= alpha; }

    // y <- alpha * (y - x)
    static void aymx(auto &y, const Scalar alpha, const auto &x) { y = alpha * (y - x); }

    // y <- alpha * x + y
    static void axpy(auto &y, const Scalar alpha, const auto &x) { y = alpha * x + y; }

    // B <- I*alpha - A
    template <typename Matrix>
    static void iama(Matrix &B, Scalar const alpha, const Matrix &A) {
        B = Matrix::Identity() * alpha - A;
    }

    // Inplace LU decomposition
    template <typename Matrix>
    static auto lu_decomposition(Matrix &A) {
        return Eigen::FullPivLU<Eigen::Ref<Matrix>>(A);
    }

    // Recalculate the LU decomposition
    template <typename Decomp, typename Matrix>
    static void update_decomposition(Decomp &decomp, Matrix &A) {
        decomp.compute(A);
    }

    // Returns `true` if matrix decomposition is non-singular
    template <typename Decomp>
    static bool invertible(const Decomp &decomp) {
        return decomp.isInvertible();
    }

    // In-place matrix solution
    template <typename Decomp, typename Vector>
    static void solve(Decomp &decomp, Vector &b) {
        b = decomp.solve(b);
    }
};


}  // namespace mathlib
}  // namespace kineticpp
