#pragma once

#include <array>
#include <type_traits>


namespace kineticpp {

template <typename T, size_t NR, size_t NC, size_t NNZ>
struct CsrMatrix
{
    using value_type = T;

    static constexpr size_t nrow = NR;
    static constexpr size_t ncol = NC;
    static constexpr size_t nnz = NNZ;

    // Structure is fixed
    const std::array<size_t, nrow+1> ridx;
    const std::array<size_t, nnz> cols;

    // Values are mutable
    std::array<value_type, nnz> vals;

    constexpr value_type operator[](size_t row, size_t col) const {
        for (size_t ii=ridx[row]; ii<ridx[row+1]; ++ii) {
            if (cols[ii] < col) continue;
            if (cols[ii] > col) return 0;
            return vals[ii];
        }
        return 0;
    }
    
    constexpr size_t rank(size_t row, size_t col) const {
        for (size_t ii=ridx[row]; ii<ridx[row+1]; ++ii) {
            if (cols[ii] < col) continue;
            if (cols[ii] > col) return nnz;
            return ii;
        }
        return nnz;
    }
};


// for (i=Start; i<End; ++i)
template <auto Start, auto End, typename B>
static constexpr void for_constexpr(B&& body)
{
    [&body]<auto... Is>(std::index_sequence<Is...>) {
        // Need std::integral_constant to make the loop index a constexpr
        (body(std::integral_constant<decltype(Start), Start+Is>()), ...);
    }(std::make_index_sequence<End-Start>{});
}

// for (i=Start; i<End; i+=Inc) 
template <auto Start, auto End, auto Inc, typename B>
static constexpr void for_constexpr(B&& body)
{
    if constexpr (Start < End) {
        // Need std::integral_constant to make the loop index a constexpr
        body(std::integral_constant<decltype(Start), Start>());
        for_constexpr<Start+Inc, End, Inc>(body);
    }
}

// for (auto x : [template_parameter_pack])
template <typename B, typename... Args>
static constexpr void foreach(B&& body, Args&&... args)
{
    (body(std::forward<Args>(args)), ...);
}
 
// for (auto x : [template_parameter_pack]) with index
template <typename B, typename... Args>
static constexpr void indexed_foreach(B&& body, Args&&... args)
{
    [&]<auto... Is>(std::index_sequence<Is...>) {
        (body(std::integral_constant<size_t, Is>(), std::forward<Args>(args)), ...);
    }(std::make_index_sequence<sizeof...(Args)>{});
}


template <auto csr, typename B>
static constexpr void for_row(B&& body) 
{
    for_constexpr<0, csr.nrow>([&](auto i) {
        body(std::integral_constant<size_t, i>());
    });
}

template <auto csr, typename B>
static constexpr void for_row_col(B&& body) 
{
    for_constexpr<0, csr.nrow>([&](auto i) {
        for_constexpr<csr.ridx[i], csr.ridx[i+1]>([&](auto ii) {
            constexpr size_t j = csr.cols[ii];
            body(std::integral_constant<size_t, i>(), std::integral_constant<size_t, j>());
        });
    });
}

template <auto csr, typename B>
static constexpr void for_nz(B&& body)
{
    for_constexpr<0, csr.nrow>([&](auto i) {
        for_constexpr<csr.ridx[i], csr.ridx[i+1]>([&](auto ii) {
            constexpr size_t j = csr.cols[ii];
            constexpr auto val = csr.vals[ii];
            body(std::integral_constant<size_t, i>(), std::integral_constant<size_t, j>(), val);
        });
    });
}

template <auto csr, size_t row, typename B>
static constexpr void for_col(B&& body)
{
    for_constexpr<csr.ridx[row], csr.ridx[row+1]>([&](auto ii) {
        constexpr size_t j = csr.cols[ii];
        body(std::integral_constant<size_t, j>());
    });
}

template <auto csr, size_t row, typename B>
static constexpr void for_nz(B&& body) 
{
    for_constexpr<csr.ridx[row], csr.ridx[row+1]>([&](auto ii) {
        constexpr size_t j = csr.cols[ii];
        constexpr auto val = csr.vals[ii];
        body(std::integral_constant<size_t, j>(), val);
    });
}


} // namespace kineticpp
