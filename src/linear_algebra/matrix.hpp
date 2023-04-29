#pragma once

#include <array>
#include <algorithm>

// FIXME namespace chem { FIXME
namespace linear_algebra {

template <typename T, size_t NR, size_t NC, size_t NNZ>
struct CsrMatrix
{
    using value_type = T;

    static constexpr size_t nnz = NNZ;
    static constexpr size_t nrow = NR;
    static constexpr size_t ncol = NC;

    const std::array<value_type, nnz> vals;
    const std::array<size_t, nnz> cols;
    const std::array<size_t, nrow+1> ridx;

    constexpr value_type operator[](size_t row, size_t col) const {
        for (size_t ii=ridx[row]; ii<ridx[row+1]; ++ii) {
            if (cols[ii] < col) continue;
            if (cols[ii] > col) return 0;
            return vals[ii];
        }
        return 0;
    }
};

} // namespace linear_algebra
//FIXME } // namespace chem
