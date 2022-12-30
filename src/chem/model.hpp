#pragma once

#include "mechanism.hpp"


namespace chem {

// template <typename LinearAlgebra>
// struct Model
// {
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
//         lhs_stoich.resize(nreact());
//         rhs_stoich.resize(nreact());
//         agg_stoich.resize(nreact());
        
//         for (size_t i=0; i<react.size(); ++i) {
//             const Reaction & rct = react[i];

//             lhs_stoich[i].resize(nspec(), 0);
//             rhs_stoich[i].resize(nspec(), 0);
//             agg_stoich[i].resize(nvar(), 0);
            
//             for (const Term & t : rct.lhs) {
//                     lhs_stoich[i][t.spc.num] += t.coef;
//                 if (_is_var_spc(t.spc)) {
//                     agg_stoich[i][t.spc.num] -= t.coef;
//                 }
//             }
//             for (const Term & t : rct.rhs) {
//                     rhs_stoich[i][t.spc.num] += t.coef;
//                 if (_is_var_spc(t.spc)) {
//                     agg_stoich[i][t.spc.num] += t.coef;
//                 }
//             }
//         }
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

//     const bool _is_var_spc(const Species & spc) {
//         return spc.num < var_spc.size();
//     }

//     const size_t nvar() const {
//         return var_spc.size();
//     }

//     const size_t nfix() const {
//         return fix_spc.size();
//     }

//     const size_t nspec() const {
//         return nvar() + nfix();
//     }

//     const size_t nreact() const {
//         return react.size();
//     }

//     vector_t calc_rates(double t) 
//     {
//         vector_t rates(mech.nreact());
//         for (size_t i=0; i<rates.size(); ++i) {
//             rates[i] = mech.react[i].rate(t);
//         }
//         return rates;
//     }

//     vector_t f(const vector_t & var, const vector_t & fix, const vector_t & rates)
//     {
//         const size_t nvar = var.size();
//         const size_t nfix = fix.size();
//         const size_t nrct = rates.size();

//         for (auto x : var) {
//             printf("%g\n", x);
//         }
//         for (auto x : fix) {
//             printf("%g\n", x);
//         }

//         vector_t vardot(nvar);
//         vector_t rate_prod(nrct);

//         printf("=====\n");
//         for (auto x : rates) {
//             printf("%g\n", x);
//         }

//         printf("=====\n");
//         for (size_t i=0; i<nrct; ++i) {
//             for (size_t j=0; j<(nvar+nfix); ++j) {
//                 printf("%4g", mech.lhs_stoich[i][j]);
//             }
//             printf("\n");
//         }
//         printf("=====\n");

//         for (size_t i=0; i<nrct; ++i) {
//             double prod = rates[i];
//             for (size_t j=0; j<nvar; ++j) {
//                 int lhs_exp = (int)(mech.lhs_stoich[i][j]);
//                 if (lhs_exp) {
//                     prod *= std::pow(var[j], lhs_exp);
//                 }
//             }
//             for (size_t j=0; j<nfix; ++j) {
//                 int lhs_exp = (int)(mech.lhs_stoich[i][nvar+j]);
//                 if (lhs_exp) {
//                     prod *= std::pow(fix[j], lhs_exp);
//                 }
//             }
//             rate_prod[i] = prod;
//         }

//         printf("-------rate_prod-------\n");
//         for (auto x : rate_prod) {
//             printf("%g\n", x);
//         }
//         printf("-------END rate_prod-------\n");

//         for (size_t j=0; j<nvar; ++j) {
//             double sum = 0;
//             for (size_t i=0; i<nrct; ++i) {    
//                 sum += mech.agg_stoich[i][j] * rate_prod[i];
//             }
//             vardot[j] = sum;
//         }

//         for (auto x : vardot) {
//             printf("%g\n", x);
//         }
//         exit(1);


//         return vardot;
//     }

//     matrix_t jac(const vector_t & var, const vector_t & fix, const vector_t & rates)
//     {
//         matrix_t jmat;
//         vector_t b(mech.lhs_stoich_rank.size());

//         const size_t nvar = var.size();
//         const size_t nfix = fix.size();
//         const size_t nrct = rates.size();

//         for (size_t i=0; i<nrct; ++i) {
//             for (size_t j=0; j<nvar; ++j) {
//                 if (mech.lhs_stoich[i][j] == 0)
//                     continue;
//                 double prod = rates[i] * mech.lhs_stoich[i][j];
//                 for (size_t k=0; k<nvar; ++k) {
//                     double m = (int)mech.lhs_stoich[i][k] - (k == j ? 1 : 0);
//                     for (size_t kk=1; kk<=m; ++kk) {
//                         prod *= var[k];
//                     }
//                 }
//                 for (size_t k=0; k<nfix; ++k) {
//                     for (size_t kk=1; kk<=(int)mech.lhs_stoich[i][nvar+k]; ++kk) {
//                         prod *= fix[k];
//                     }
//                 }
//                 b[mech.lhs_stoich_rank[i][j]-1] = prod;
//             }
//         }

//         jmat.resize(nvar);
//         for (size_t i=0; i<nvar; ++i) {
//             jmat[i].resize(nvar);
//             for (size_t j=0; j<nvar; ++j) {
//                 if (!mech.jac_struct[i][j])
//                     continue;
//                 double sum = 0;
//                 for (size_t k=0; k<nrct; ++k) {
//                     if (mech.agg_stoich[k][i] && mech.lhs_stoich_rank[k][j]) {
//                         sum += mech.agg_stoich[k][i] * b[mech.lhs_stoich_rank[k][j]-1];
//                     }
//                 }
//                 jmat[i][j] = sum;
//             }
//         }

//         return jmat;
//     }

// };


} // namespace chem