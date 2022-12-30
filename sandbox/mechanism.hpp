#pragma once
#include <array>
#include <utility>
#include <type_traits>
#include "atoms.hpp"
#include "util.hpp"


// Chemical species atomic composition
template <typename... Atoms>
struct Species
{
    static constexpr size_t mass = (Atoms::mass + ...);
};


template <typename... Spec>
struct SpeciesGroup
{
    static constexpr size_t size = sizeof...(Spec);
};


template <auto C, typename S>
struct Term
{
    static constexpr double coef = C;
    using spec = S;
};


template <typename... Terms>
struct Expression
{
    template <typename>
    static constexpr bool involves() {
        return false;
    }
    
    template <typename Species, typename Car, typename... Cdr>
    static constexpr bool involves() {
        if constexpr (is_instance<Car, Term>()) {
            return std::is_same_v<Species, Car::S> || involves<Species, Cdr...>();
        } else {
            return std::is_same_v<Species, Car> || involves<Species, Cdr...>();
        }
    }

};

template <typename LHS, typename RHS, auto Rate>
struct Reaction
{
    template <typename Species>
    static constexpr bool lhs_involves() {
        return false;
    }
};

template <
    typename T,
    typename VarSpecies, 
    typename FixSpecies, 
    typename... Reactions>
struct Mechanism
{
    static constexpr size_t nvar = VarSpecies::size;
    static constexpr size_t nfix = FixSpecies::size;
    static constexpr size_t nspec = nvar + nfix;
    static constexpr size_t nreact = sizeof...(Reactions);

    // Stoichiometric matrix: [nvar+nfix, nreact]
    using stoich_mat_t = std::array<std::array<T, nreact>, nspec>;

    static constexpr stoich_mat_t _lhs_stoichiometry(VarSpecies&& var, FixSpecies&& fix, Reactions&&... react) {
        stoich_mat_t lhs_stoich;
        // for_constexpr(size_t spc=0; spc<NVAR; ++spc) {
        //     for_constexpr([&](auto const & rct) {
        //         lhs_stoich[1][1] = 1;
        //     }, react...);
        // }
        return lhs_stoich;
    }
    static constexpr stoich_mat_t lhs_stoich = _lhs_stoichiometry(VarSpecies(), FixSpecies(), Reactions()...);
};

