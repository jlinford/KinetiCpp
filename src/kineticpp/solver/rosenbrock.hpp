// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <cmath>
#include <limits>

#include "status.hpp"


namespace kineticpp {
namespace solver {


template <typename P, typename LA, typename VarConc, typename FixConc, typename Jacobian>
struct Rosenbrock {

    static constexpr double double_max = std::numeric_limits<double>::max();

    struct Args {
        ErrorCode status;          // Integration status
        size_t maxstep = 100000;   // Maximum number of solver steps
        size_t maxdecomp = 5;      // Maximum number of attempted Jacobian matrix decompositions
        bool scrub = true;         // Clamp near-zero values to zero
        double abstol = 1.0;       // Absolute tolerance
        double reltol = 1e-3;      // Relative tolerance
        double h = 0;              // Integration timestep
        double hmin = 0;           // Minimum timestep size
        double hmax = double_max;  // Maximum timestep size
        double hlim = 1e-5;        // Lower bound on step size
        double hfac = 0.5;         // Step size reduction factor in case of failed Jacobian matrix decomp
        double facmin = 0.2;       // Lower bound on new step calculation factor
        double facmax = 6;         // Upper bound on new step calculation factor
        double facrej = 0.1;       // Step size back up rejection factor
        double facsafe = 0.9;      // Safety factor for new step size calculation
        size_t nstep = 0;          // Number of integration steps
        size_t nfun = 0;           // Number of ODE function calls
        size_t njac = 0;           // Number of ODE Jacobian calls
        size_t ndecomp = 0;        // Number of matrix decompositions
        size_t nsolve = 0;         // Number of matrix solutions
    };

    static ErrorCode integrate(auto fun, auto jac, VarConc &var, const FixConc &fix, const double t0, const double tend,
                               Args &args) {

        // Initialize linear solver
        auto solver = typename LA::Solver();

        // Allocate ODE function vectors
        std::array<VarConc, P::S> K;  // Stage vectors
        VarConc f0;                   // ODE function f(t)
        VarConc dfdt;                 // Finite difference approximation of df/dt
        VarConc fs;                   // ODE function evaulation for stage
        VarConc udot;                 // New solution
        VarConc uerr;                 // New solution error

        // Allocate compressed Jacobian nonzero vectors
        Jacobian j0;     // Jacobian fJ(t)
        Jacobian hgimj;  // Step matrix: hgimj = 1/(h*gamma)*I - fJ(t)

        // Current integration time
        double t = t0;

        // Initial integration timestep
        args.h = std::min(std::max(args.hmin, std::max(args.hmin, args.hlim)), std::max(args.hmax, tend - t0));

        // Time integration
        while (t < tend) {
            if (args.nstep++ > args.maxstep) {
                return (args.status = ErrorCode::iterations);
            }
            args.h = std::min(args.h, std::abs(tend - t));

            fun(f0, var, fix, t);
            ++args.nfun;

            jac(j0, var, fix, t);
            ++args.njac;

            // Finite difference approximation of df/dt
            constexpr double eps = 10 * std::numeric_limits<double>::epsilon();
            constexpr double delmin = 1e-6;
            const double tdel = std::sqrt(eps) * std::max(delmin, std::abs(t));
            fun(dfdt, var, fix, t + tdel);
            ++args.nfun;
            LA::aymx(dfdt, 1.0 / tdel, f0);

            // Step calculation
            bool reject = false;
            while (true) {
                // Construct step matrix: hgimj = 1/(h*gamma)*I - fJ(t)
                LA::iama(hgimj, 1.0 / (args.h * P::Gamma[0]), j0);

                // Calculate LU decomposition of step matrix
                auto decomp_succcess = solver.decompose(hgimj);
                ++args.ndecomp;

                // If decomposition fails, reduce step size and retry
                if (!decomp_succcess) {
                    double hbar = args.h * args.hfac;
                    for (size_t ndecomp = 1; ndecomp < args.maxdecomp; ++ndecomp) {
                        LA::iama(hgimj, 1.0 / (hbar * P::Gamma[0]), j0);
                        auto decomp_succcess = solver.decompose(hgimj);
                        ++args.ndecomp;
                        if (decomp_succcess) {
                            break;
                        } else if (ndecomp == args.maxdecomp - 1) {
                            return (args.status = ErrorCode::decomposition);
                        }
                        hbar *= args.hfac;
                    }
                    args.h = hbar;
                }

                // Calculate stage vectors
                for (size_t stage = 0; stage < P::S; ++stage) {
                    VarConc &sK = K[stage];
                    if (stage == 0) {
                        LA::copy(fs, f0);
                    } else {
                        if (P::EvalF[stage]) {
                            VarConc ubar;
                            LA::copy(ubar, var);
                            for (size_t i = 0; i < stage; ++i) {
                                double alpha = P::A[stage * (stage - 1) / 2 + i];
                                LA::axpy(ubar, alpha, K[i]);
                            }
                            double tau = t + args.h * P::Alpha[stage];
                            fun(fs, ubar, fix, tau);
                            ++args.nfun;
                        }
                    }
                    LA::copy(sK, fs);
                    for (size_t i = 0; i < stage; ++i) {
                        double Ch = P::C[stage * (stage - 1) / 2 + i] / args.h;
                        LA::axpy(sK, Ch, K[i]);
                    }
                    if (P::Gamma[stage]) {
                        double hGamma = args.h * P::Gamma[stage];
                        LA::axpy(sK, hGamma, dfdt);
                    }
                    solver.solve(sK);
                    ++args.nsolve;
                }  // Stage calculation

                // New solution
                LA::copy(udot, var);
                for (size_t i = 0; i < P::S; ++i) {
                    LA::axpy(udot, P::M[i], K[i]);
                }

                // Estimate error
                LA::zero(uerr);
                for (size_t i = 0; i < P::S; ++i) {
                    LA::axpy(uerr, P::E[i], K[i]);
                }
                double err = 0;
                for (auto i = 0; i < LA::size(var); ++i) {
                    double umax = std::max(std::abs(var[i]), std::abs(udot[i]));
                    double scale = args.abstol + args.reltol * umax;
                    err += (uerr[i] * uerr[i]) / (scale * scale);
                }
                err = std::sqrt(err / LA::size(var));

                // New step size
                double hdot =
                    args.h * std::min(args.facmax, std::max(args.facmin, args.facsafe / std::pow(err, 1.0 / P::ELO)));
                if ((err <= 1.0) || (args.h < args.hmin)) {
                    // Accept step and update solution
                    LA::copy(var, udot, args.scrub);
                    t += args.h;
                    // Set step size for next iteration
                    hdot = std::max(args.hmin, std::min(hdot, std::max(args.hmax, tend - t0)));
                    args.h = reject ? std::min(hdot, args.h) : hdot;
                    // Exit step calcuation loop
                    break;
                }
                if (reject) {
                    // Consecutive rejections
                    // Back off step size by rejection factor and repeat step
                    // calculation
                    args.h *= args.facrej;
                } else {
                    // First rejection
                    // Retry step calculation with the new step size before
                    // backing off step size
                    reject = true;
                    args.h = hdot;
                }
            }  // while (true)
        }      // while (t <= tend)

        // Integration successful
        return (args.status = ErrorCode::success);
    }  // Integrate

};  // class Rosenbrock


namespace rosenbrock_parameters {

struct Ros2 {
    // Precision of integration parameters
    using T = double;

    // Estimator of local order
    static constexpr T ELO = 2.0;

    // Number of stages
    static constexpr int S = 2;

    // True: Evaluate f()
    // False: Reuse f() from previous stage
    static constexpr std::array<bool, S> EvalF {true, true};

    // Lower triangular coefficient matrix
    // A[i,j] = A[i*(i-1)/2 + j]
    static constexpr std::array<T, 1> A {0.5857864376};

    // Lower triangular coefficient matrix
    // C[i,j] = C[i*(i-1)/2 + j]
    static constexpr std::array<T, 1> C {-1.1715728753};

    // New step coefficients
    static constexpr std::array<T, S> M {
        0.8786796564,
        0.2928932188,
    };

    // Error estimator coefficients
    static constexpr std::array<T, 2> E {0.2928932188, 0.2928932188};

    // Y_i ~= Y(T+H*Alpha[i])
    static constexpr std::array<T, 2> Alpha {0.0, 1.0};

    // Gamma_i = sum(Gamma[i,:])
    static constexpr std::array<T, 2> Gamma {1.7071067812, -1.7071067812};
};


struct Ros4 {
    // Precision of integration parameters
    using T = double;

    // Estimator of local order
    static constexpr T ELO = 4.0;

    // Number of stages
    static constexpr int S = 4;

    // True: Evaluate f()
    // False: Reuse f() from previous stage
    static constexpr std::array<bool, S> EvalF {true, true, true, false};

    // Lower triangular coefficient matrix
    // A[i,j] = A[i*(i-1)/2 + j]
    static constexpr std::array<T, 6> A {
        2.0, 1.867943637803922, 0.2344449711399156, 1.867943637803922, 0.2344449711399156, 0};

    // Lower triangular coefficient matrix
    // C[i,j] = C[i*(i-1)/2 + j]
    static constexpr std::array<T, 6> C {-7.137615036412310, 2.580708087951457,   0.6515950076447975,
                                         -2.137148994382534, -0.3214669691237626, -0.6949742501781779};

    // New step coefficients
    static constexpr std::array<T, S> M {2.255570073418735, 0.2870493262186792, 0.4353179431840180, 1.093502252409163};

    // Error estimator coefficients
    static constexpr std::array<T, S> E {-0.2815431932141155, -0.07276199124938920, -0.1082196201495311,
                                         -1.093502252409163};

    // Y_i ~= Y(T+H*Alpha[i])
    static constexpr std::array<T, S> Alpha {
        0.0,
        1.145640000000000,
        0.6552168638155900,
        0.6552168638155900,
    };

    // Gamma_i = sum(Gamma[i,:])
    static constexpr std::array<T, S> Gamma {0.5728200000000000, -1.769193891319233, 0.7592633437920482,
                                             -0.1049021087100450};
};

}  // namespace rosenbrock_parameters


template <typename... Ts>
using Ros2 = Rosenbrock<rosenbrock_parameters::Ros2, Ts...>;

template <typename... Ts>
using Ros4 = Rosenbrock<rosenbrock_parameters::Ros4, Ts...>;


}  // namespace solver
}  // namespace kineticpp