#pragma once

namespace kineticpp {

struct Expression {};

template <typename L, typename R>
struct BinaryOperator : public Expression
{ 
    const L lhs;
    const R rhs;
    constexpr BinaryOperator(L& lhs, R& rhs) :
        lhs(lhs), rhs(rhs)
    {}
};

template <typename L, typename R>
struct Sum : public BinaryOperator<L,R>
{ 
    constexpr Sum(L& lhs, R& rhs) : 
        BinaryOperator<L,R>(lhs, rhs)
    {}
};

template <typename L, typename R>
struct Product : public BinaryOperator<L,R>
{ 
    constexpr Product(L& lhs, R& rhs) : 
        BinaryOperator<L,R>(lhs, rhs)
    {}
};

template <typename L, typename R>
struct Exponential : public BinaryOperator<L,R>
{ 
    constexpr Exponential(L& lhs, R& rhs) : 
        BinaryOperator<L,R>(lhs, rhs)
    {}
};

template <typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

template <typename T>
concept Term = (std::is_arithmetic_v<T> || std::is_enum_v<T> || std::is_base_of_v<Expression, T>);

template <Term L, Term R>
static constexpr auto operator+(L lhs, R rhs) {
    return Sum {lhs, rhs};
}

template <Term L, Term R>
static constexpr auto operator*(L lhs, R rhs) {
    return Product {lhs, rhs};
}

template <Term L, Term R>
static constexpr auto operator^(L lhs, R rhs) {
    return Exponential {lhs, rhs};
}

} // namespace kineticpp