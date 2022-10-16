#pragma once

template <std::size_t Start, std::size_t End, typename Body>
constexpr void for_constexpr(Body && body)
{
    [&body]<size_t... Is>(std::index_sequence<Is...>) {
        // Need std::integral_constant to make the loop index a constexpr
        // (body(Start + Is), ...);
        (body(std::integral_constant<decltype(Start), Start+Is>()), ...);
    }(std::make_index_sequence<End-Start>{});
}


template <typename Body, typename... Args>
constexpr void for_constexpr(Body && body, Args && ... args)
{
    (body(std::forward<Args>(args)), ...);
}