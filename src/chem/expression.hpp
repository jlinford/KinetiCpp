#pragma once

namespace chem {

template <typename L, typename R>
struct Sum
{
    const L lhs;
    const R rhs;
};

template <typename L, typename R>
struct Product
{
    const L lhs;
    const R rhs;
};

template <typename L, typename R>
struct Exponential
{
    const L lhs;
    const R rhs;
};

template <typename L, typename R>
static constexpr auto operator+(L lhs, R rhs) {
    return Sum {lhs, rhs};
}

template <typename L, typename R>
static constexpr auto operator*(L lhs, R rhs) {
    return Product {lhs, rhs};
}

template <typename L, typename R>
static constexpr auto operator^(L lhs, R rhs) {
    return Exponential {lhs, rhs};
}

template <typename T>
concept Arithmetic = std::is_arithmetic_v<T>;


}