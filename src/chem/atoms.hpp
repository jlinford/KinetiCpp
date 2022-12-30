#pragma once


namespace chem {

struct Atom 
{
    unsigned int number;
};

// Zero mass balance: ignore atomic mass balance
static constexpr Atom ZMB {0};

// Atom definitions for mass balance calculation
static constexpr Atom H {1};
static constexpr Atom He {2};
static constexpr Atom Li {3};
static constexpr Atom Be {4};
static constexpr Atom B {5};
static constexpr Atom C {6};
static constexpr Atom N {7};
static constexpr Atom O {8};
static constexpr Atom F {9};
static constexpr Atom Ne {10};

} // namespace Chem
