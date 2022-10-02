#pragma once
#include <cstddef>
#include <array>

template <typename T, size_t NNZ, size_t M, size_t N>
struct CSR
{
    using elem_t = T;
    using vals_t = std::array<T, NNZ>; 
    using cidx_t = std::array<const size_t, NNZ>;
    using ridx_t = std::array<const size_t, M+1>;

    constexpr CSR(vals_t const & _vals, cidx_t const & _cidx, ridx_t const & _ridx) :
            vals(_vals), cidx(_cidx), ridx(_ridx)
    { }

    constexpr size_t row_start(size_t i) const {
        return ridx[i];
    }

    constexpr size_t row_end(size_t i) const {
        return ridx[i+1];
    }

    vals_t vals;
    const cidx_t cidx;
    const ridx_t ridx;
};