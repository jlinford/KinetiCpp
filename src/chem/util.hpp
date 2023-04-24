#pragma once

#include <utility>


// template <auto Start, auto End, typename B>
// static constexpr void for_constexpr(B&& body)
// {
//     [&body]<auto... Is>(std::index_sequence<Is...>) {
//         // Need std::integral_constant to make the loop index a constexpr
//         (body(std::integral_constant<decltype(Start), Start+Is>()), ...);
//     }(std::make_index_sequence<End-Start>{});
// }

template <auto Start, auto End, typename B>
static constexpr void for_constexpr(B&& body)
{
    [&body]<auto... Is>(std::index_sequence<Is...>) {
        // Need std::integral_constant to make the loop index a constexpr
        (body(std::integral_constant<decltype(Start), Start+Is>()), ...);
    }(std::make_index_sequence<End-Start>{});
}

// template <auto Start, auto End, auto Inc, typename B>
// static constexpr void for_constexpr(B&& body)
// {
//     if constexpr (Start < End) {
//         body(std::integral_constant<decltype(Start), Start>());
//         for_constexpr<Start+Inc, End, Inc>(body);
//     }
// }

template <auto Start, auto End, auto Inc, typename B>
static constexpr void for_constexpr(B&& body)
{
    if constexpr (Start < End) {
        body(std::integral_constant<decltype(Start), Start>());
        for_constexpr<Start+Inc, End, Inc>(body);
    }
}

// template <typename B, typename... Args>
// static constexpr void for_constexpr(B&& body, Args&&... args)
// {
//     (body(std::forward<Args>(args)), ...);
// }

template <typename B, typename... Args>
static constexpr void for_constexpr(B&& body, Args&&... args)
{
    (body(std::forward<Args>(args)), ...);
}

// template <class B, class T>
// static constexpr void for_constexpr(B&& body, T&& tuple)
// {
//     for_constexpr<size_t(0), std::tuple_size_v<std::decay_t<T>>, size_t(1)>([&](auto i) {
//         body(std::get<i.value>(tuple));
//     });
// }

template <class B, class T>
static constexpr void for_constexpr(B&& body, T&& tuple)
{
    for_constexpr<size_t(0), std::tuple_size_v<std::decay_t<T>>, size_t(1)>([&](auto i) {
        body(std::get<i.value>(tuple));
    });
}












template <auto Start, auto End, auto Inc, class F>
constexpr void constexpr_for(F&& f)
{
    if constexpr (Start < End)
    {
        f(std::integral_constant<decltype(Start), Start>());
        constexpr_for<Start + Inc, End, Inc>(f);
    }
}

template <class F, class Tuple>
constexpr void constexpr_for_tuple(F&& f, Tuple&& tuple)
{
    constexpr size_t cnt = std::tuple_size_v<std::decay_t<Tuple>>;

    constexpr_for<size_t(0), cnt, size_t(1)>([&](auto i) {
        f(std::get<i.value>(tuple));
    });
}