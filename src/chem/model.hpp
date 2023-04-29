#pragma once

#include <chem/mechanism.hpp>
#include <solver/rosenbrock.hpp>


namespace chem {

template <typename M, typename LA>
class Model
{
public:

    template <size_t N>
    using vector_t = typename LA:: template Vector<N>;

    template <size_t Rows, size_t Cols>
    using matrix_t = typename LA:: template Matrix<Rows, Cols>;

    template <typename Solver, typename... Args>
    static void solve(auto & conc, double t0, double tend, Args... arg) {
        //solver::RosenbrockImpl<M::nspc, Solver, LA>::integrate(fun, jac, conc, t0, tend, arg...);
    }

// private:

    template <typename T> requires std::is_floating_point_v<T>
    static constexpr bool is_nonzero(T elem) {
        return std::abs(elem) > (10*std::numeric_limits<T>::epsilon());
    }
    
    static void fun(vector_t<M::nspc>& du, const vector_t<M::nspc>& u, const double t) {
        constexpr auto lhs = M::lhs_stoich();
        constexpr auto agg = M::agg_stoich();
        auto rates = M::rates(t);
        vector_t<M::nrct> rate_prod;
        for_constexpr<0, lhs.nrow>([&](auto i) {
            double prod = rates[i];
            for_constexpr<lhs.ridx[i], lhs.ridx[i+1]>([&](auto ii) {
                size_t j = lhs.cols[ii];
                double lhs_exp = lhs.vals[ii];
                prod *= std::pow(u[j], lhs_exp);
            });
            rate_prod[i] = prod;
        });
        for_constexpr<0, M::nvar>([&](auto j) {
            double sum = 0;
            for_constexpr<0, M::nrct>([&](auto i) {
                constexpr double agg_coef = agg[i,j];
                if constexpr (is_nonzero(agg_coef)) {
                    sum += agg_coef * rate_prod[i];
                }
            });
            du[j] = sum;
        });
        // for_constexpr<0, M::nvar>([&](auto j) {
        //     double sum = 0;
        //     for_constexpr<0, M::nrct>([&](auto i) {
        //         constexpr double agg_coef = agg[i,j];
        //         if constexpr (is_nonzero(agg_coef)) {
        //             sum += agg_coef * rate_prod[i];
        //         }
        //     });
        //     du[j] = sum;
        // });
        for_constexpr<0, M::nfix>([&](auto j) {
            du[M::nvar+j] = 0;
        });
    }

    template <auto lhs_stoich, auto agg_stoich>
    static void jac(matrix_t<M::nspc, M::nspc>& J, const vector_t<M::nspc>& u, const double t) {
        constexpr auto lhs = std::experimental::mdspan(lhs_stoich.data(), M::nrct, M::nspc);
        constexpr auto agg = std::experimental::mdspan(agg_stoich.data(), M::nrct, M::nvar);

        auto rates = M::rates(t);
        // matrix_t<M::nrct, M::nvar> B;
        std::array<double, M::nrct*M::nvar> B;
        B.fill(0);

        for_constexpr<0, M::nrct>([&](auto i) {
            for_constexpr<1, 1+M::nvar>([&](auto j) {
                if constexpr (is_nonzero(lhs[i,j-1])) {
                    double p = rates[i];
                    // p *= prod(u[1:j-1].^lhs_stoich[1:j-1,i])
                    for_constexpr<0, j-1>([&](auto k) {
                        constexpr double lhs_exp = lhs[i,k];
                        if constexpr (is_nonzero(lhs_exp)) {
                            p *= std::pow(u[k], lhs_exp);
                        }
                    });
                    // p *= lhs_stoich[j,i] * u[j]^(lhs_stoich[j,i]-1)
                    double lhs_exp = lhs[i,j-1];
                    p *= lhs_exp * std::pow(u[j-1], lhs_exp-1);
                    // p *= prod(u[j+1:nspec].^lhs_stoich[j+1:nspec,i])
                    for_constexpr<j.value, M::nspc>([&](auto k) {
                        constexpr double lhs_exp = lhs[i,k];
                        if constexpr (is_nonzero(lhs_exp)) {
                            p *= std::pow(u[k], lhs_exp);
                        }
                    });
                    B[i*M::nvar+j-1] = p;
                }
            });
        });

        // std::array<double, M::nvar*M::nvar> Jtmp;
        for (size_t k=0; k<M::nrct; ++k) {
            for (size_t i=0; i<M::nvar; ++i) {
                for (size_t j=0; j<M::nvar; ++j) {
                    if (i == j || is_nonzero(agg[k,i]*lhs[k,j])) {
                        // J[i,j] = sum(agg_stoich[:,i] .* B[:,j])
                        double sum = 0;
                        for (size_t ii=0; ii<M::nrct; ++ii) {
                            //sum += agg[ii,i] * B[ii,j];
                            sum += agg[ii,i] * B[ii*M::nvar+j];
                        }
                        // Jtmp[i*M::nvar+j] = sum;
                        J(i,j) = sum;
                    }
                }
            }
        }
    }


//     using Vector = typename LinearAlgebra:: template Vector<N>;
//     using Matrix = typename LinearAlgebra:: template Matrix<N, N>;
//     using LUDecomp = typename LinearAlgebra:: template LUDecomp<Matrix>;


//     // Stoichiometric matrices
//     // Row indices are reactions
//     // Column indices are variable species
//     matrix_t<double> lhs_stoich;
//     matrix_t<double> rhs_stoich;
//     matrix_t<double> agg_stoich;
//     // Used in computing intermediate values of the Jacobian
//     matrix_t<size_t> lhs_stoich_rank;

//     // Species Jacobian sparsity structure
//     matrix_t<bool> jac_struct;


//     Mechanism(
//         const species_list_t && _var_spc, 
//         const species_list_t && _fix_spc,
//         const reaction_list_t && _react) :
//             var_spc(_var_spc),
//             fix_spc(_fix_spc),
//             react(_react)
//         {
//             _index_species();
//             _construct_stoichiometric_matrix();
//             _construct_jacobian_sparsity();
//         }

//     void _index_species() 
//     {
//         size_t num = 0;
//         for (Species & spc : var_spc) {
//             spc.num = num;
//             ++num;
//         }
//         for (Species & spc : fix_spc) {
//             spc.num = num;
//             ++num;
//         }
//     }
    
//     void _construct_stoichiometric_matrix()
//     {

//         size_t nz = 0;
//         lhs_stoich_rank.resize(nreact());
//         for (size_t i=0; i<nreact(); ++i) {
//             lhs_stoich_rank[i].resize(nspec());
//             for (size_t j=0; j<nspec(); ++j) {
//                 if (lhs_stoich[i][j]) {
//                     ++nz;
//                     lhs_stoich_rank[i][j] = nz;
//                 }
//             }
//         }
//     }

//     void _construct_jacobian_sparsity()
//     {
//         jac_struct.resize(nvar());
//         for (size_t i=0; i<nvar(); ++i) {
//             jac_struct[i].resize(nvar());
//             jac_struct[i][i] = true;
//         }
//         for (size_t i=0; i<nvar(); ++i) {
//             for (size_t j=0; j<nvar(); ++j) {
//                 for (size_t k=0; k<nreact(); ++k) {
//                     if (lhs_stoich[j][k] && agg_stoich[i][k]) {
//                         jac_struct[i][j] = true;
//                     }
//                 }
//             }
//         }
//     }

}; // Model

} // namespace chem