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

static constexpr bool is_zero(auto elem) { return !is_nonzero(elem); }

template <size_t... Is>
constexpr auto reverse_index_sequence(std::index_sequence<Is...> const &)
    -> decltype(std::index_sequence<sizeof...(Is) - 1U - Is...> {});

template <size_t N>
using make_reverse_index_sequence = decltype(reverse_index_sequence(std::make_index_sequence<N> {}));

// for (i=Start; i<End; ++i) or for(i=End-1; i>=Start; --i)
// std::integral_constant makes the loop index a constexpr
template <auto Start, auto End>
static constexpr void for_constexpr(auto &&body) {
    if constexpr (Start < End) {
        [&body]<auto... Is>(std::index_sequence<Is...>) {
            (body(std::integral_constant<size_t, Start + Is>()), ...);
        }(std::make_index_sequence<End - Start> {});
    } else if constexpr (End < Start) {
        [&body]<auto... Is>(std::index_sequence<Is...>) {
            (body(std::integral_constant<size_t, End + Is>()), ...);
        }(make_reverse_index_sequence<Start - End> {});
    }
}

// for (i=Start; i<End; i+=Inc)
template <auto Start, auto End, auto Inc>
static constexpr void for_constexpr(auto &&body) {
    if constexpr (Start < End) {
        body(std::integral_constant<size_t, Start>());
        for_constexpr<Start + Inc, End, Inc>(body);
    }
}

// for (auto x : [template_parameter_pack])
template <typename... Args>
static constexpr void foreach(auto &&body, Args &&...args) {
    (body(std::forward<Args>(args)), ...);
}

// for (auto x : [template_parameter_pack]) with index
template <typename... Args>
static constexpr void indexed_foreach(auto &&body, Args &&...args) {
    [&]<auto... Is>(std::index_sequence<Is...>) {
        (body(std::integral_constant<size_t, Is>(), std::forward<Args>(args)), ...);
    }(std::make_index_sequence<sizeof...(Args)> {});
}

}  // namespace kineticpp