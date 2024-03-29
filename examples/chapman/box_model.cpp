//
// A small Chapman-like mechanism:
// https://glossary.ametsoc.org/wiki/Chapman_mechanism
//
// SPDX-License-Identifier: Apache-2.0
// Copyright 2023, John Linford <john@redhpc.com>
//

#include <cmath>
#include <iostream>

#include <kineticpp/math/stdcpp.hpp>
#include <kineticpp/model.hpp>
#include <kineticpp/solver/rosenbrock.hpp>

#include "chapman.hpp"


// Program entry point
int main(int argc, char **argv) {

    // Define a model of a chemical mechanism.
    // Use the Chapman-like mechanism with a Ros4 implicit time stepping
    // integrator implemented with Eigen dense matrices
    using Model = kineticpp::Model<double, photo_chem::Chapman, kineticpp::solver::Ros4, kineticpp::math::StdCpp>;

    // Initial concentrations.
    // Species order is specified at mechanism definition
    Model::VarConc var {
        9.906E+01,  // O1D
        6.624E+08,  // O1
        5.326E+11,  // O3
        8.725E+08,  // NO
        2.240E+08   // NO2
    };
    Model::FixConc fix {
        8.120E+16,  // M
        1.697E+16   // O2
    };

    double h = 0;             // solver timestep on exit
    double t0 = 0;            // integration start time
    double tend = 24 * 3600;  // integration end time
    double dt = 600;          // time between output

    // output header
    std::cout << "step;t;h;O1D;O1;O3;NO;NO2;M;O2" << std::endl;

    // Time integration
    // Use callback to report concentrations
    auto errcode = Model::solve(var, fix, t0, tend, dt, [](auto &step, auto &t, auto &var, auto &fix, auto &args) {
        std::cout << step << ";" << t << ";" << args.h << ";";
        for (auto &x : var) {
            std::cout << x << ";";
        }
        for (auto &x : fix) {
            std::cout << x << ";";
        }
        std::cout << std::endl;
    });

    // How'd it go?
    std::cout << kineticpp::solver::explain(errcode) << std::endl;

    return (errcode != kineticpp::solver::ErrorCode::success) ? EXIT_FAILURE : EXIT_SUCCESS;
}