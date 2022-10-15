#pragma once
#include <cstddef>

template <size_t Start, size_t End, typename Body>
constexpr void for_constexpr(Body && body)
{
    [&body]<size_t... Is>(std::index_sequence<Is...>) {
        (body(Start + Is), ...);
    }(std::make_index_sequence<End-Start>{});
}
