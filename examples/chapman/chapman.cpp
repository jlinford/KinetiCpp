//
// A small Chapman-like mechanism: 
// https://glossary.ametsoc.org/wiki/Chapman_mechanism
//
// SPDX-License-Identifier: Apache-2.0
// Copyright 2023, John Linford <john@redhpc.com>
//

#include <iostream>
#include <cmath>
#include <kineticpp/model.hpp>
#include <kineticpp/mathlib/eigen.hpp>


// A simple sunlight intensity model
double sunlight(double seconds) 
{
    const double sunrise = 4.5 * 3600;
    const double sunset  = 19.5 * 3600;
    if (seconds < sunrise || seconds > sunset) {
        return 0;
    }
    
    const double tmp = std::abs((2.0*seconds-sunrise-sunset)/(sunset-sunrise));
    return (1.0 + std::cos(M_PI*(tmp*tmp))) / 2.0;
}


int main() 
{   
    using namespace kineticpp;
    using namespace kineticpp::atom;

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
        O2 >= 2*O1          || [](double t) { return 2.643e-10 * std::pow(sunlight(t), 3); },
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

    using Model = Model<SmallChapman, solver::Ros4, mathlib::EigenDense<double>>;
    
    Model::Solution conc {
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

    double h = 0;
    double t0 = 0;
    double tend = 24*3600;
    double dt = 600;

    std::cout << "step;t;h;O1D;O1;O3;NO;NO2;M;O2";
    auto errcode = Model::solve(conc, h, t0, tend, dt, [](auto& step, auto& t, auto& h, auto& u) {
        std::cout << step << ";" << t << ";" << h << ";";
        for (auto& x : u) {
            std::cout << x << ";";
        }
        std::cout << std::endl;
    });

    std::cout << solver::explain(errcode) << std::endl;

    return (errcode != solver::ErrorCode::success) ? EXIT_FAILURE : EXIT_SUCCESS;
}