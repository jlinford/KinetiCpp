// Copyright 2023, John Linford <john@redhpc.com>
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <cmath>
#include <limits>

#include "status.hpp"


namespace kineticpp {
namespace solver {


template <typename P, size_t N, typename LA>
class Rosenbrock {
public:
    using LinearAlgebra = LA;
    using Vector = typename LA::template Vector<N>;
    using Matrix = typename LA::template Matrix<N, N>;
    using LUDecomp = typename LA::template LUDecomp<Matrix>;

    struct Parameters {
        double abstol = 1.0;
        double reltol = 1e-3;
        size_t maxstep = 100000;
        size_t maxdecomp = 5;
        double hmin = 0;
        double hmax = std::numeric_limits<double>::max();
        double hlim = 1e-5;
        double hfac = 0.5;
        double facmin = 0.2;
        double facmax = 6;
        double facrej = 0.1;
        double facsafe = 0.9;
    };

    static ErrorCode integrate(auto fun, auto jac, Vector &u, double &h, const double t0, const double tend,
                               const Parameters &args) {

        constexpr double eps = 10 * std::numeric_limits<double>::epsilon();
        constexpr double delmin = 1e-6;

        // Integration step count
        size_t nstep = 0;

        // Current integration time
        double t = t0;

        // Integration timestep
        h = std::min(std::max(args.hmin, std::max(args.hmin, args.hlim)), std::max(args.hmax, tend - t0));

        // Stage vectors
        std::array<Vector, P::S> K;

        // Time integration
        while (t < tend) {
            if (nstep++ > args.maxstep) {
                return ErrorCode::iterations;
            }
            h = std::min(h, std::abs(tend - t));

            Vector f0;
            fun(f0, u, t);

            Matrix j0;
            jac(j0, u, t);

            // Finite difference approximation of df/dt
            const double tdel = std::sqrt(eps) * std::max(delmin, std::abs(t));
            const double tdel_inv = 1.0 / tdel;
            Vector dfdt;
            fun(dfdt, u, t + tdel);
            LA::aymx(dfdt, tdel_inv, f0);

            // Step calculation
            bool reject = false;
            while (true) {
                // Construct step matrix: hgimj = 1/(h*gamma)*I - fJ(t)
                Matrix hgimj;
                LA::iama(hgimj, 1.0 / (h * P::Gamma[0]), j0);

                // Calculate LU decomposition of step matrix
                LUDecomp hgimj_decomp = LA::lu_decomposition(hgimj);

                // If decomposition is not invertable reduce step size and try
                // again
                if (!LA::invertible(hgimj_decomp)) {
                    double hbar = h * args.hfac;
                    for (size_t ndecomp = 1; ndecomp < args.maxdecomp; ++ndecomp) {
                        LA::iama(hgimj, 1.0 / (hbar * P::Gamma[0]), j0);
                        LA::update_decomposition(hgimj_decomp, hgimj);
                        if (LA::invertible(hgimj_decomp)) {
                            // Decomp successful
                            break;
                        } else if (ndecomp == args.maxdecomp - 1) {
                            return ErrorCode::decomposition;
                        }
                        hbar *= args.hfac;
                    }
                    h = hbar;
                }

                // Calculate stages
                for (size_t stage = 0; stage < P::S; ++stage) {
                    Vector &sK = K[stage];
                    Vector fs;
                    if (stage == 0) {
                        LA::copy(fs, f0);
                    } else {
                        if (P::EvalF[stage]) {
                            Vector ubar;
                            LA::copy(ubar, u);
                            for (size_t i = 0; i < stage; ++i) {
                                double alpha = P::A[stage * (stage - 1) / 2 + i];
                                LA::axpy(ubar, alpha, K[i]);
                            }
                            double tau = t + h * P::Alpha[stage];
                            fun(fs, ubar, tau);
                        }
                    }
                    LA::copy(sK, fs);
                    for (size_t i = 0; i < stage; ++i) {
                        double Ch = P::C[stage * (stage - 1) / 2 + i] / h;
                        LA::axpy(sK, Ch, K[i]);
                    }
                    if (P::Gamma[stage]) {
                        double hGamma = h * P::Gamma[stage];
                        LA::axpy(sK, hGamma, dfdt);
                    }
                    LA::solve(hgimj_decomp, sK);
                }  // Stage calculation

                // New solution
                Vector udot;
                LA::copy(udot, u);
                for (size_t i = 0; i < P::S; ++i) {
                    LA::axpy(udot, P::M[i], K[i]);
                }

                // Estimate error
                Vector uerr;
                LA::zero(uerr);
                for (size_t i = 0; i < P::S; ++i) {
                    LA::axpy(uerr, P::E[i], K[i]);
                }
                double err = 0;
                for (decltype(u.size()) i = 0; i < u.size(); ++i) {
                    double umax = std::max(std::abs(u[i]), std::abs(udot[i]));
                    double scale = args.abstol + args.reltol * umax;
                    err += (uerr[i] * uerr[i]) / (scale * scale);
                }
                err = std::sqrt(err / u.size());

                // New step size
                double hdot =
                    h * std::min(args.facmax, std::max(args.facmin, args.facsafe / std::pow(err, 1.0 / P::ELO)));
                if ((err <= 1.0) || (h < args.hmin)) {
                    // Accept step and update solution
                    LA::copy(u, udot);
                    t += h;
                    // Set step size for next iteration
                    hdot = std::max(args.hmin, std::min(hdot, std::max(args.hmax, tend - t0)));
                    h = reject ? std::min(hdot, h) : hdot;
                    // exit loop
                    break;
                }
                if (reject) {
                    // Consecutive rejections
                    // Back off step size by rejection factor and repeat step
                    // calculation
                    h *= args.facrej;
                } else {
                    // First rejection
                    // Retry step calculation with the new step size before
                    // backing off step size
                    reject = true;
                    h = hdot;
                }
            }  // while (true)
        }      // while (t <= tend)

        // Integration successful
        return ErrorCode::success;
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


template <size_t N, typename LA>
using Ros2 = Rosenbrock<rosenbrock_parameters::Ros2, N, LA>;

template <size_t N, typename LA>
using Ros4 = Rosenbrock<rosenbrock_parameters::Ros4, N, LA>;


}  // namespace solver
}  // namespace kineticpp