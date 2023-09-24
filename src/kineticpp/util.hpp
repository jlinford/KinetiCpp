// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <limits>
#include <type_traits>

namespace kineticpp {

template <typename T>
requires std::is_floating_point_v<T>
static constexpr bool is_nonzero(T elem) {
    return std::abs(elem) > (10 * std::numeric_limits<T>::epsilon());
}

// for (i=Start; i<End; ++i)
template <auto Start, auto End, typename B>
static constexpr void for_constexpr(B &&body) {
    [&body]<auto... Is>(std::index_sequence<Is...>) {
        // Need std::integral_constant to make the loop index a constexpr
        (body(std::integral_constant<decltype(Start), Start + Is>()), ...);
    }(std::make_index_sequence<End - Start> {});
}

// for (i=Start; i<End; i+=Inc)
template <auto Start, auto End, auto Inc, typename B>
static constexpr void for_constexpr(B &&body) {
    if constexpr (Start < End) {
        // Need std::integral_constant to make the loop index a constexpr
        body(std::integral_constant<decltype(Start), Start>());
        for_constexpr<Start + Inc, End, Inc>(body);
    }
}

// for (auto x : [template_parameter_pack])
template <typename B, typename... Args>
static constexpr void foreach(B &&body, Args &&...args) {
    (body(std::forward<Args>(args)), ...);
}

// for (auto x : [template_parameter_pack]) with index
template <typename B, typename... Args>
static constexpr void indexed_foreach(B &&body, Args &&...args) {
    [&]<auto... Is>(std::index_sequence<Is...>) {
        (body(std::integral_constant<size_t, Is>(), std::forward<Args>(args)), ...);
    }(std::make_index_sequence<sizeof...(Args)> {});
}

}  // namespace kineticpp