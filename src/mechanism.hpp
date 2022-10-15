#pragma once
#include <array>
#include "atoms.hpp"


// Chemical species atomic composition
template <typename... Atoms>
struct Species
{

};

template <typename... Spec>
struct SpeciesGroup
{
    static constexpr size_t SIZE = sizeof...(Spec);
};

template <auto N, typename Spec>
struct Term
{

};

template <typename... Terms>
struct Expression
{};

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
        
    };
};

