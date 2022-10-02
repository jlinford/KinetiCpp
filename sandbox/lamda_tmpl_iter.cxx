#include <cstddef>
#include <utility>
#include <array>



template <size_t Start, size_t End, typename Body>
constexpr void for_constexpr(Body && body)
{
    [&body]<size_t... Is>(std::index_sequence<Is...>) {
        (body(std::integral_constant<decltype(Start), Start+Is>()), ...);
    }(std::make_index_sequence<End-Start>{});
}

template <size_t Start, size_t End, signed long Inc, typename Body>
constexpr void for_constexpr(Body && body)
{
    static_assert(Inc != 0);

    if constexpr (Inc > 0) {
        if constexpr (Start < End) {
            body(std::integral_constant<decltype(Start), Start>());
            for_constexpr<Start+Inc, End, Inc>(body);
        }
    } else {
        if constexpr (Start >= End) {
            body(Start);
            for_constexpr<Start+Inc, End, Inc>(body);
        }
    }
}


int main(int argc, char ** argv)
{
    constexpr std::array<double, 5> nz{1, 2, 3, 4, 5};

    for_constexpr<1, 3>([&](auto i) {
        constexpr double value = nz[i];
    });

    return 0;
}