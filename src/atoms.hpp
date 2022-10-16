#pragma once

namespace atoms {

using atomic_number_t = unsigned int;

template <atomic_number_t N>
struct Atom {
    static constexpr atomic_number_t num = N;
};

// Ignore atomic mass balance
using IGNORE = Atom<0>;

// Atom definitions for mass balance calculation
using H = Atom<1>;
using He = Atom<2>;
using Li = Atom<3>;
using Be = Atom<4>;
using B = Atom<5>;
using C = Atom<6>;
using N = Atom<7>;
using O = Atom<8>;
using F = Atom<9>;
using Ne = Atom<10>;

} // namespace atoms
