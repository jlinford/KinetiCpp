// #include <iostream>
#include <algorithm>
#include <array>

//#include "csr.hpp"
template <typename T,
          size_t NNZ, 
          size_t M, 
          size_t N,
          std::array<size_t,NNZ> CIDX,
          std::array<size_t,M+1> RIDX>
struct CSR
{
    using elem_t = T;
    using nz_t = std::array<T, NNZ>; 

    CSR(nz_t const & _nz) : nz(_nz) {}

    constexpr static size_t nnz = NNZ;
    constexpr static size_t rows = M;
    constexpr static size_t cols = N;
    constexpr static auto cidx = CIDX;
    constexpr static auto ridx = RIDX;

    const nz_t nz;
};


//#include "for_constexpr.hpp"
template <size_t Start, size_t End, typename Body>
constexpr void for_constexpr(Body && body)
{
    [&body]<size_t... Is>(std::index_sequence<Is...>) {
        // Need std::integral_constant to make the loop index a constexpr
        // (body(Start + Is), ...);
        (body(std::integral_constant<decltype(Start), Start+Is>()), ...);
    }(std::make_index_sequence<End-Start>{});
}


// template <typename T>
// void print_container(T const & v)
// {
//     for (auto x : v) {
//         std::cout << x << ", ";
//     }
//     std::cout << std::endl;
// }


template <typename T, typename M>
void sp_gemv_constexpr(
    std::array<T, M::rows> & y, 
    M & A,
    std::array<T, M::cols> x, 
    T const alpha=1.0,
    T const beta=1.0)
{
    // y := alpha*A*x + beta*y,
    for_constexpr<0, M::rows>([&](auto const i) {
        constexpr size_t rstart = A.ridx[i];
        constexpr size_t rend = A.ridx[i+1];
        T sum = 0;
        for_constexpr<rstart, rend>([&](auto const j) {
            constexpr size_t col = A.cidx[j];
            const T value = A.nz[j];
            sum += value * x[col];
        });
        y[i] = beta*y[i] + alpha*sum;
    });
}


int main(int argc, char ** argv)
{
    // MxN matrix with 4 nonzeros
    //     A          x
    // [ 1 0 0 ]   [  0 ]   [  0 ]
    // [ 0 2 0 ]   [ 10 ]   [ 20 ]
    // [ 0 0 3 ] x [ 20 ] = [ 60 ]
    // [ 0 4 0 ]   [ 30 ]   [ 40 ]
    constexpr size_t M = 4;
    constexpr size_t N = 3;
    constexpr size_t NNZ = 4;

    constexpr std::array<size_t, NNZ> cidx {0, 1, 2, 1};
    constexpr std::array<size_t, M+1> ridx {0, 1, 2, 3, 4};

    std::array<double, NNZ> nz;
    // for (size_t i=0; i<nz.size(); ++i) {
    //     nz[i] = i+1;
    // }
    // print_container(nz);

    std::array<double, N> x;
    for(size_t i=0; i<x.size(); ++i) {
        x[i] = i*10;
    }
    // print_container(x);

    std::array<double, M> y{};

    CSR<double, NNZ, M, N, cidx, ridx> A(nz);

    sp_gemv_constexpr<double>(y, A, x);
    // print_container(y);

    double sum = 0;
    for (auto x : y) {
        sum += x;
    };

    return sum;
}