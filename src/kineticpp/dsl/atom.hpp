// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cmath>
#include <type_traits>

#include <kineticpp/dsl/expression.hpp>


namespace kineticpp::dsl {

// clang-format off
enum AtomicSymbol : size_t
{
    // Atomic number 0 has zero mass, i.e. ignore mass balance
    ZERO_MASS = 0,
    // Atomic numbers
     H,                                                                                                                         He,
    Li, Be,                                                                                                  B,  C,  N,  O,  F, Ne,
    Na, Mg,                                                                                                 Al, Si,  P,  S, Cl, Ar,
     K, Ca, Sc,                                                         Ti,  V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
    Rb, Sr,  Y,                                                         Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te,  I, Xe,
    Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta,  W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra, Ac, Th, Pa,  U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og
};
// clang-format on

static constexpr double AtomicMass[] {
    0,       1.007,   4.002,  6.941,   9.012,  10.811,  12.011,  14.007,  15.999,  18.998,  20.18,   22.99,
    24.305,  26.982,  28.086, 30.974,  32.065, 35.453,  39.948,  39.098,  40.078,  44.956,  47.867,  50.942,
    51.996,  54.938,  55.845, 58.933,  58.693, 63.546,  65.38,   69.723,  72.64,   74.922,  78.96,   79.904,
    83.798,  85.468,  87.62,  88.906,  91.224, 92.906,  95.96,   98,      101.07,  102.906, 106.42,  107.868,
    112.411, 114.818, 118.71, 121.76,  127.6,  126.904, 131.293, 132.905, 137.327, 138.905, 140.116, 140.908,
    144.242, 145,     150.36, 151.964, 157.25, 158.925, 162.5,   164.93,  167.259, 168.934, 173.054, 174.967,
    178.49,  180.948, 183.84, 186.207, 190.23, 192.217, 195.084, 196.967, 200.59,  204.383, 207.2,   208.98,
    210,     210,     222,    223,     226,    227,     232.038, 231.036, 238.029, 237,     244,     243,
    247,     247,     251,    252,     257,    258,     259,     262,     261,     262,     266,     264,
    267,     268,     271,    272,     285,    284,     289,     288,     292,     295,     294};


static constexpr double atomic_mass(AtomicSymbol atm) { return AtomicMass[static_cast<size_t>(atm)]; }

template <typename T>
static constexpr bool is_atom(const T &) {
    return std::is_same_v<std::decay_t<T>, AtomicSymbol>;
}

template <typename T>
requires std::is_arithmetic_v<T>
static constexpr double atomic_mass(T atm) {
    return atm;
}

template <typename L, typename R>
static constexpr double atomic_mass(Sum<L, R> term) {
    return atomic_mass(term.lhs) + atomic_mass(term.rhs);
}

template <typename L, typename R>
static constexpr double atomic_mass(Product<L, R> term) {
    return atomic_mass(term.lhs) * atomic_mass(term.rhs);
}

template <typename L, typename R>
static constexpr double atomic_mass(Exponential<L, R> term) {
    return std::pow(atomic_mass(term.lhs), atomic_mass(term.rhs));
}

}  // namespace kineticpp::dsl
