// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <array>

#include <kineticpp/util.hpp>


namespace kineticpp::math {

template <typename T, typename JacStruct, typename JacLUStruct>
struct StdCpp {

    template <size_t N>
    using Vector = std::array<T, N>;

    using Scalar = T;
    using Jacobian = Vector<JacStruct::nnz>;

    template <size_t N>
    static constexpr size_t size(Vector<N> &y) {
        return y.size();
    }

    // y <- 0
    template <size_t N>
    static void zero(Vector<N> &y) {
        y.fill(0);
    }

    // y <- x
    template <size_t N>
    static void copy(Vector<N> &y, const Vector<N> &x, bool scrub = false) {
        if (scrub) {
            std::transform(x.cbegin(), x.cend(), y.begin(), [](const auto &x_i) {
                return (util::is_nonzero(x_i) ? x_i : 0);
            });
        } else {
            std::copy(x.cbegin(), x.cend(), y.begin());
        }
    }

    // y <- alpha * y
    template <size_t N>
    static void scale(Vector<N> &y, const Scalar alpha) {
        std::transform(y.cbegin(), y.cend(), y.begin(y), [alpha](const auto &y_i) {
            return alpha * y_i;
        });
    }

    // y <- alpha * (y - x)
    template <size_t N>
    static void aymx(Vector<N> &y, const Scalar alpha, const Vector<N> &x) {
        std::transform(x.cbegin(), x.cend(), y.cbegin(), y.begin(), [alpha](const auto &x_i, const auto &y_i) {
            return alpha * (y_i - x_i);
        });
    }

    // y <- alpha * x + y
    template <size_t N>
    static void axpy(Vector<N> &y, const Scalar alpha, const Vector<N> &x) {
        std::transform(x.cbegin(), x.cend(), y.cbegin(), y.begin(), [alpha](const auto &x_i, const auto &y_i) {
            return alpha * x_i + y_i;
        });
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

    class Solver {
        Vector<JacLUStruct::nnz> jac_lu;

    public:
        bool decompose(Jacobian &J) {
            for (size_t i = 0; i < JacStruct::nrow; ++i) {
                if (util::is_zero(J[JacStruct::diag(i)])) {
                    return false;
                }
            }
            Vector<JacLUStruct::ncol> row;
            JacLUStruct::for_ridx([&](auto i) {
                JacLUStruct::for_cidx_rank_in_row(i, [&](auto j, auto rank) {
                    if constexpr (JacStruct::value(i, j)) {
                        row[j] = J[JacStruct::rank(i, j)];
                    } else {
                        row[j] = 0;
                    }
                });
                JacLUStruct::for_cidx_rank_below_diag(i, [&](auto j, auto rank) {
                    double a = -row[j] / jac_lu[JacLUStruct::diag(j)];
                    row[j] = -a;
                    JacLUStruct::for_cidx_rank_above_diag(j, [&](auto k, auto rank) {
                        row[k] += a * jac_lu[rank];
                    });
                });
                JacLUStruct::for_cidx_rank_in_row(i, [&](auto j, auto rank) {
                    jac_lu[rank] = row[j];
                });
            });
            return true;
        }

        bool solve(Vector<JacLUStruct::nrow> &x) {
            JacLUStruct::for_ridx([&](auto i) {
                if constexpr (JacLUStruct::ridx(i) < JacLUStruct::diag(i)) {
                    double sum = x[i];
                    JacLUStruct::for_cidx_rank_below_diag(i, [&](auto j, auto rank) {
                        sum -= jac_lu[rank] * x[j];
                    });
                    x[i] = sum;
                }
            });
            JacLUStruct::for_ridx_reversed([&](auto i) {
                double sum = x[i];
                JacLUStruct::for_cidx_rank_above_diag(i, [&](auto j, auto rank) {
                    sum -= jac_lu[rank] * x[j];
                });
                x[i] = sum / jac_lu[JacLUStruct::diag(i)];
            });
            return true;
        }
    };
};

}  // namespace kineticpp::math
