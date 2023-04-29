#include <cstdio>
#include <cmath>
#include <linear_algebra/eigen.hpp>
#include <chem/model.hpp>

// A simple sunlight model
double sunlight(double t) 
{
    const double sunrise = 4.5 * 3600;
    const double sunset  = 19.5 * 3600;
    if (t < sunrise || t > sunset) return 0;
    
    double tmp = std::abs((2.0*t-sunrise-sunset)/(sunset-sunrise));
    return (1.0 + std::cos(M_PI*(tmp*tmp))) / 2.0;
}


int main() 
{   
    // A translation of the "Small Strato" mechanism from KPP (Sandu et al.)
    // <R1>  O2   + hv = 2O		    : (2.643E-10) * SUN*SUN*SUN;
    // <R2>  O    + O2 = O3		    : (8.018E-17);
    // <R3>  O3   + hv = O   + O2 	: (6.120E-04) * SUN;
    // <R4>  O    + O3 = 2O2		: (1.576E-15);
    // <R5>  O3   + hv = O1D + O2	: (1.070E-03) * SUN*SUN;
    // <R6>  O1D  + M  = O   + M	: (7.110E-11);
    // <R7>  O1D  + O3 = 2O2 		: (1.200E-10);
    // <R8>  NO   + O3 = NO2 + O2 	: (6.062E-15);
    // <R9>  NO2  + O  = NO  + O2	: (1.069E-11);
    // <R10> NO2  + hv = NO  + O	: (1.289E-02) * SUN;

    using namespace chem;
    using namespace chem::atom;

    // Declare chemical species names
    enum S {
        O1D=100, O1, O3, NO, NO2, M, O2
    };

    using SmallChapman = Mechanism<
        // Variable species concentrations change according to the law of mass action kinetics
        // Order determines species numbering
        VariableSpecies {
            O1D || O,
            O1  || O,
            O3  || O*3,
            NO  || N+O,
            NO2 || N+O*2
        },
        // Fixed species concentrations are determined by physical factors
        // Order determines species numbering
        FixedSpecies {
            M  || N*2 + O*2,
            O2 || O*2
        },
        // Reactions and rates
        // Order determines equation numbering
        O2 >= 2*O1          || [](double t) { return 2.643e-10 * std::pow(sunlight(t), 3); }
        ,
        O1 + O2 >= O3       || 8.018e-17,
        O3 >= O1 + O2       || [](double t) { return 6.12e-04 * sunlight(t); },
        O1 + O3 >= 2*O2     || 1.576e-15,
        O3 >= O1D + O2      || [](double t) { return 1.07e-03 * std::pow(sunlight(t), 2); },
        O1D + M >= O1 + M   || 7.11e-11,
        O1D + O3 >= 2*O2    || 1.2e-10,
        NO + O3 >= NO2 + O2 || 6.062e-15,
        NO2 + O1 >= NO + O2 || 1.069e-11,
        NO2 >= NO + O1      || [](double t) { return 1.289e-02 * sunlight(t); }
    >;

    using LinearAlgebra = linear_algebra::EigenLib<double>;
    using Model = chem::Model<SmallChapman, LinearAlgebra>;
    
    Model::vector_t<SmallChapman::nspc> conc {
        // Initial concentrations of variable species
        9.906E+01,
        6.624E+08,
        5.326E+11,
        8.725E+08,
        2.240E+08,
        // Initial concentrations of fixed species
        8.120E+16,
        1.697E+16
    };

    // auto mat = std::experimental::mdspan(Model::agg_stoich.data(), SmallChapman::nrct, SmallChapman::nvar);
    // for (size_t i=0; i<mat.extent(0); ++i) {
    //     for (size_t j=0; j<mat.extent(1); ++j) {
    //         if (mat[i,j]) {
    //             std::cout << mat[i,j] << " ";
    //         } else {
    //             std::cout << "." << " ";
    //         }
    //     }
    //     std::cout << std::endl;
    // }

    double tt = 12*3600;
    Model::vector_t<SmallChapman::nspc> du;
    //Model::fun<SmallChapman::lhs_stoich(), SmallChapman::agg_stoich()>(du, conc, tt);
    Model::fun(du, conc, tt);
    for (auto x: du) {
        printf("%g\n", x);
    }

    // double tt = 12*3600;
    // Model::matrix_t<SmallChapman::nspc, SmallChapman::nspc> J;
    // Model::jac<SmallChapman::lhs_stoich(), SmallChapman::agg_stoich()> (J, conc, tt);
    // for (int i=0; i<J.rows(); ++i) {
    //     for (int j=0; j<J.cols(); ++j) {
    //         printf("%16g  ", J(i,j));
    //         //std::cout << J(i,j) << "  ";
    //     }
    //     //std::cout << std::endl;
    //     printf("\n");
    // }

    // Model::solve<solver::Ros4>(conc, 0, 24*3600);

    return 0;
}