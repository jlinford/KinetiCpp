#pragma once
#include <array>
#include "atoms.hpp"
#include "util.hpp"


// Chemical species atomic composition
template <typename... Atoms>
struct Species
{
    static constexpr size_t weight = 0;
};

template <typename... Spec>
struct SpeciesGroup
{
    static constexpr size_t SIZE = sizeof...(Spec);
};

template <auto X, typename S>
struct Term
{
    static constexpr auto coef = X;
    using Spec = S;
};

template <typename... Terms>
struct Expression
{

};

template <typename LHS, typename RHS, auto Rate>
struct Reaction
{};

template <
    typename T,
    typename VarSpecies, 
    typename FixSpecies, 
    typename... Reactions>
struct Mechanism
{
    static constexpr size_t NVAR = VarSpecies::SIZE;
    static constexpr size_t NFIX = FixSpecies::SIZE;
    static constexpr size_t NSPEC = NVAR + NFIX;
    static constexpr size_t NREACT = sizeof...(Reactions);

    // Stoichiometric matrix: [NSPEC, NREACT]
    using stoich_mat_t = std::array<std::array<T, NREACT>, NSPEC>;
    static constexpr stoich_mat_t stoich = []() {
        stoich_mat_t stoich;
        // for_constexpr([](auto const & rct) {

        // }, react...);
        return stoich;
    };
};

