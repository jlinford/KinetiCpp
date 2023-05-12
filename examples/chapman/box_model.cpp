//
// A small Chapman-like mechanism:
// https://glossary.ametsoc.org/wiki/Chapman_mechanism
//
// SPDX-License-Identifier: Apache-2.0
// Copyright 2023, John Linford <john@redhpc.com>
//

#include <cmath>
#include <iostream>
#include <kineticpp/mathlib/eigen.hpp>
#include <kineticpp/model.hpp>
#include <kineticpp/solver/rosenbrock.hpp>

#include "chapman.hpp"


// Program entry point
int main(int argc, char **argv) {

    // Define a model of a chemical mechanism.
    // Use the Chapman-like mechanism with a Ros4 implicit time stepping
    // integrator implemented with Eigen dense matrices
    using Model =
        kineticpp::Model<photo_chem::Chapman, kineticpp::solver::Ros4, kineticpp::mathlib::EigenDense<double>>;

    // Initial concentrations.
    // Species order is specified at mechanism definition
    Model::Solution conc {
        9.906E+01,  // O1D
        6.624E+08,  // O1
        5.326E+11,  // O3
        8.725E+08,  // NO
        2.240E+08,  // NO2
        8.120E+16,  // M
        1.697E+16   // O2
    };

    double h = 0;             // solver timestep on exit
    double t0 = 0;            // integration start time
    double tend = 24 * 3600;  // integration end time
    double dt = 600;          // time between output

    // output header
    std::cout << "step;t;h;O1D;O1;O3;NO;NO2;M;O2";

    // Time integration
    // Use callback to report concentrations
    auto errcode = Model::solve(
        [](auto &step, auto &t, auto &h, auto &u) {
            std::cout << step << ";" << t << ";" << h << ";";
            for (auto &x : u) {
                std::cout << x << ";";
            }
            std::cout << std::endl;
        },
        conc, h, t0, tend, dt);

    // How'd it go?
    std::cout << kineticpp::solver::explain(errcode) << std::endl;

    return (errcode != kineticpp::solver::ErrorCode::success) ? EXIT_FAILURE : EXIT_SUCCESS;
}