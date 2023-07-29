//
// Performance benchmarks
//
// SPDX-License-Identifier: Apache-2.0
// Copyright 2023, John Linford <john@redhpc.com>
//

#include <chrono>
#include <cmath>
#include <iostream>
#include <kineticpp/mathlib/eigen.hpp>
#include <kineticpp/model.hpp>
#include <kineticpp/solver/rosenbrock.hpp>

#include "chapman.hpp"


// Program entry point
int main(int argc, char **argv) {

    constexpr size_t iterations = 10000;

    auto benchmark = [&](auto &label, auto &&body) {
        auto tstart = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < iterations; ++i) {
            body();
        }
        auto tend = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::micro> duration = tend - tstart;
        std::cout << label << ": " << iterations << " calls in " << duration << ".  " << (duration / iterations)
                  << " avg per call" << std::endl;
    };

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

    double t = 12 * 3600;

    Model::Solver::Vector du;
    benchmark("Function", [&]() { Model::fun(du, conc, t); });

    Model::Solver::Matrix J;
    benchmark("Jacobian", [&]() { Model::jac(J, conc, t); });

    benchmark("Time Integration", [&]() {
        double h = 0;
        Model::solve(conc, h, 0, 24 * 3600);
    });

    return EXIT_SUCCESS;
}