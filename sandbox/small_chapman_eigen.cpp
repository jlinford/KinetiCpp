#include <cstdio>
#include <cmath>
#include <chem/model.hpp>
#include <linear_algebra/eigen.hpp>



// A simple sunlight model
double sunlight(double t) 
{
    const double sunrise = 4.5 * 3600;
    const double sunset  = 19.5 * 3600;
    if (t < sunrise || t > sunset) {
        return 0;
    }
    double tmp = std::abs((2.0*t-sunrise-sunset)/(sunset-sunrise));
    return (1.0 + std::cos(M_PI*(tmp*tmp))) / 2.0;
}


int main(int argc, char ** argv) 
{   
    using Atom = chem::Atom;
    const Atom _O = Atom::O;
    const Atom _N = Atom::N;

    // Declare chemical species
    using Species = chem::Species;
    Species O {{_O}};
    Species O1D {{_O}};
    Species O3 {{_O, _O, _O}};
    Species NO {{_N, _O}};
    Species NO2 {{_N, _O, _O}};
    Species M {{_N, _N, _O, _O}};
    Species O2 {{_O, _O}};

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
    using Mechanism = chem::Mechanism<5, 2, 10>;
    Mechanism small_chapman {
        // Variable species concentrations change according to the law of mass action kinetics
        // Order determines species numbering
        {O1D, O, O3, NO, NO2},
        // Fixed species concentrations are determined by physical factors
        // Order determines species numbering
        {M, O2},
        // Reactions and rates
        // Order determines equation numbering
        {
            { // O2 + hv -> 2*O 
                {O2},
                {{2.0,O}},
                [](double t) {
                    double sun = sunlight(t);
                    return 2.643e-10 * sun*sun*sun;
                }
            },
            { // O + O2 -> O3
                {O, O2},
                {O3},
                [](double) {
                    return 8.018e-17;
                }
            },
            { // O3 + hv -> O + O2
                {O3},
                {O, O2},
                [](double t) {
                    double sun = sunlight(t);
                    return 6.12e-04 * sun;
                }
            },
            { // O + O3 -> 2*O2
                {O, O3},
                {{2.0,O2}},
                [](double) {
                    return 1.576e-15;
                }
            },
            { // O3 + hv -> O1D + O2
                {O3},
                {O1D, O2},
                [](double t) {
                    double sun = sunlight(t);
                    return 1.07e-03 * sun*sun;
                }
            },
            { // O1D + M -> O + M
                {O1D,M},
                {O,M},
                [](double) {
                    return 7.11e-11;
                }
            },
            { // O1D + O3 -> 2*O2
                {O1D, O3},
                {{2.0,O2}},
                [](double) {
                    return 1.2e-10;
                }
            },
            { // NO + O3 -> NO2 + O2
                {NO, O3},
                {NO2, O2},
                [](double) {
                    return 6.062e-15;
                }
            },
            { // NO2 + O -> NO + O2
                {NO2, O},
                {NO, O2},
                [](double) {
                    return 1.069e-11;
                }
            },
            { // NO2 + hv -> NO + O
                {NO2},
                {NO, O},
                [](double t) {
                    double sun = sunlight(t);
                    return 1.289e-02 * sun;
                }
            }
        }
    };

    using LinearAlgebra = linear_algebra::EigenLib<double>;
    using Model = chem::Model<Mechanism, LinearAlgebra>;
    
    Model small_model(small_chapman); 
    
    Model::Vector conc {
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

    small_model.solve<solver::Ros4>(conc, 0, 24*3600);

    return 0;
}