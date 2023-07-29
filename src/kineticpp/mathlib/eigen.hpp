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

    static constexpr size_t size(auto &y) {
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
    static void iama(auto &B, Scalar const alpha, const auto &A) {
        using Matrix = std::remove_reference_t<decltype(A)>;
        B = Matrix::Identity() * alpha - A;
    }

    // Inplace LU decomposition
    static auto lu_decomposition(auto &A) {
        using Matrix = std::remove_reference_t<decltype(A)>;
        return Eigen::FullPivLU<Eigen::Ref<Matrix>>(A);
    }

    // Recalculate the LU decomposition
    static void update_decomposition(auto &decomp, auto &A) {
        decomp.compute(A);
    }

    // Returns `true` if matrix decomposition is non-singular
    static bool invertible(const auto &decomp) {
        return decomp.isInvertible();
    }

    // In-place matrix solution
    static void solve(auto &decomp, auto &b) {
        b = decomp.solve(b);
    }
};


}  // namespace mathlib
}  // namespace kineticpp
