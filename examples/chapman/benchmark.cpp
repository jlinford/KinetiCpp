//
// Performance benchmarks
//
// SPDX-License-Identifier: Apache-2.0
// Copyright 2023, John Linford <john@redhpc.com>
//

#include <array>
#include <chrono>
#include <cmath>
#include <iostream>

#include <kineticpp/mathlib/eigen.hpp>
#include <kineticpp/mathlib/stdcpp.hpp>
#include <kineticpp/model.hpp>
#include <kineticpp/solver/rosenbrock.hpp>

#include "chapman.hpp"


template <typename Model>
void benchmark_model() {

    constexpr size_t iterations = 50000;

    auto benchmark = [&](auto &label, auto &&body) {
        std::cout << label << ": " << std::flush;
        auto tstart = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < iterations; ++i) {
            body();
        }
        auto tend = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::micro> duration = tend - tstart;
        std::cout << iterations << " calls in " << duration << ".  " << (duration / iterations) << " avg per call"
                  << std::endl;
    };

    // Initial concentrations.
    // Species order is specified at mechanism definition
    typename Model::VarConc var {
        9.906E+01,  // O1D
        6.624E+08,  // O1
        5.326E+11,  // O3
        8.725E+08,  // NO
        2.240E+08   // NO2
    };
    typename Model::FixConc fix {
        8.120E+16,  // M
        1.697E+16   // O2
    };

    double t = 12 * 3600;

    typename Model::VarConc du;
    benchmark("Function", [&]() {
        Model::fun(du, var, fix, t);
    });

    typename Model::Jacobian J;
    benchmark("Jacobian", [&]() {
        Model::jac(J, var, fix, t);
    });

    benchmark("Integrate", [&]() {
        Model::solve(var, fix, 0, 24 * 3600);
    });
}


int main(int argc, char **argv) {

    std::cout << "Eigen" << std::endl;
    benchmark_model<
        kineticpp::Model<double, photo_chem::Chapman, kineticpp::solver::Ros4, kineticpp::mathlib::Eigen>>();

    std::cout << "StdCpp" << std::endl;
    benchmark_model<
        kineticpp::Model<double, photo_chem::Chapman, kineticpp::solver::Ros4, kineticpp::mathlib::StdCpp>>();

    return EXIT_SUCCESS;
}
