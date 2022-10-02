#include <iostream>
#include <array>
#include <algorithm>
#include <cstdlib>
#include <utility>


// https://artificial-mind.net/blog/2020/10/31/constexpr-for
// template <auto Start, auto End, auto Inc, class F>
// constexpr void constexpr_for(F && f)
// {
//     if constexpr (Start < End) {
//         f(std::integral_constant<decltype(Start), Start>());
//         constexpr_for<Start+Inc, End, Inc>(f);
//     }
// }

template <auto Start, auto End, auto Inc, template<auto> class F>
constexpr void constexpr_for()
{
    if constexpr (Start < End) {
        F<Start>();
    }
}


template <size_t istart, size_t iend, typename T, size_t NNZ, size_t NIDX>
constexpr T row_sum(const std::array<T,NNZ> & x, const std::array<size_t,NIDX> & idx)
{
    T sum = 0;
    // constexpr_for<istart, iend, 1>([&](auto i) {
    //     const size_t j = idx[i];
    //     sum += x[j];
    // });

    [&]<size_t i...>(std::index_sequence<i...>) {
        const size_t j = idx[i];
        sum += x[j];
    }
    (std::make_index_sequence<iend-istart+1>{});

    // constexpr_for<istart, iend, 1, body>();

    return sum;
}


int main(int argc, char ** argv)
{
    size_t const NNZ = 8;
    size_t const NIDX = 4;

    std::array<size_t, NIDX> idx {0, 6, 4, 2};
    
    // std::array<int, NNZ> x {10, 20, 30, 40, 50, 60, 70, 80};
    
    std::array<int, NNZ> x;
    // std::generate(x.begin(), x.end(), rand);
    
    return row_sum<1,4>(x, idx);
}