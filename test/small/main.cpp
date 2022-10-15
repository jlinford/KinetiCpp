#include <iostream>

#include "mechanism.hpp"
#include "model.hpp"


int main(void)
{
    using atoms::IGNORE;
    using _O = atoms::O;
    using _N = atoms::N;

    // Declare chemical species
    using O = Species<_O>;
    using O1D = Species<_O>;
    using O3 = Species<_O, _O, _O>;
    using NO = Species<_N, _O>;
    using NO2 = Species<_N, _O, _O>;
    using M = Species<_N, _N, _O, _O>;
    using O2 = Species<_O, _O>;

    // A fake species to represent radiative sunlight
    using hv = Species<IGNORE>;
    
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
    using SmallChapman = Mechanism<
        // Precision of stoichiometric data structures
        double,
        // Variable species concentrations change according to the law of mass action kinetics
        SpeciesGroup<O, O1D, O3, NO, NO2>,
        // Fixed species concentrations are determined by physical factors
        SpeciesGroup<M, O2>,
        // Reactions and rates
        Reaction<
            Expression<O2, hv>,
            Expression<Term<2.0,O>>,
            [](double sun) {
                return 2.643e-10 * sun*sun*sun;
            }
        >,
        Reaction<
            Expression<O, O2>,
            Expression<O3>,
            8.018e-17
        >,
        Reaction<
            Expression<O3, hv>,
            Expression<O, O2>,
            [](double sun) {
                return 6.120E-04 * sun;
            }
        >,
        Reaction<
            Expression<O, O3>,
            Expression<Term<2.0,O2>>,
            1.576e-15
        >,
        Reaction<
            Expression<O3, hv>,
            Expression<O1D,O2>,
            [](double sun) {
                return 1.070e-03 * sun*sun;
            }
        >,
        Reaction<
            Expression<O1D,M>,
            Expression<O,M>,
            7.110E-11
        >,
        Reaction<
            Expression<O1D,O3>,
            Expression<Term<2.0,O2>>,
            1.200E-10
        >,
        Reaction<
            Expression<NO, O3>,
            Expression<NO2, O2>,
            6.062E-15
        >,
        Reaction<
            Expression<NO2, O>,
            Expression<NO, O2>,
            1.069E-11
        >,
        Reaction<
            Expression<NO2, hv>,
            Expression<NO, O>,
            [](double sun) {
                return 1.289e-02 * sun;
            }
        >
    >; // SmallChapman

    using SmallModel = Model<
        // Precision of concentration data
        double, 
        // Integration method
        Ros2, 
        // Chemical mechanism
        SmallChapman>;

    // Integration time grid
    double t0 = 12*3600;
    double tend = 24*3600;
    double tdel = 0.25*3600;

    // Integration tolerances
    double abstol = 1.0;
    double reltol = 1e-3;

    // Initial concentrations of variable species
    SmallModel::var_t var_conc {
        9.906E+01,
        6.624E+08,
        5.326E+11,
        8.725E+08,
        2.240E+08
    };

    // Initial concentrations of fixed species
    SmallModel::fix_t fix_conc {
        8.120E+16,
        1.697E+16
    };

    // Time integration
    SmallModel().integrate(
        t0, tend, tdel, 
        var_conc, fix_conc,
        abstol, reltol);

    return 0;
}