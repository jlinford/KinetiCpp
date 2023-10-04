// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>
#include <type_traits>

#include <kineticpp/dsl/atom.hpp>


namespace kineticpp::dsl {

template <typename T>
concept SpeciesId = std::is_enum_v<T>;

template <SpeciesId ID>
struct Species {
    using type = Species<ID>;
    using id_type = ID;

    const id_type id;
    const double mass;

    constexpr operator id_type() const { return static_cast<id_type>(id); }
};

template <typename L, typename R>
struct Equation {
    const L lhs;
    const R rhs;
};

template <typename L, typename R, typename Rate>
struct Reaction {
    const Equation<L, R> eqn;
    const Rate rate;
};

template <typename ID, Term T>
static constexpr auto operator||(ID lhs, T rhs) {
    return Species<ID> {lhs, atomic_mass(rhs)};
}

// template <Term L, Term R>
template <typename L, typename R>
static constexpr auto operator>=(L lhs, R rhs) {
    return Equation<L, R> {lhs, rhs};
}

template <typename L, typename R, typename Rate>
static constexpr auto operator||(Equation<L, R> lhs, Rate rhs) {
    return Reaction<L, R, Rate> {lhs, rhs};
}

}  // namespace kineticpp::dsl