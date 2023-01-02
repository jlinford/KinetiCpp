#pragma once

#include <cmath>
#include <array>
#include <limits>


namespace solver {


struct Ros2
{
    // Precision of integration parameters
    using T = double;

    // Estimator of local order
    static constexpr T ELO = 2.0;
    
    // Number of stages
    static constexpr int S = 2;

    // True: Evaluate f()
    // False: Reuse f() from previous stage
    static constexpr std::array<bool, S> EvalF {
        true, 
        true
    };

    // Lower triangular coefficient matrix
    // A[i,j] = A[i*(i-1)/2 + j]
    static constexpr std::array<T, 1> A {
        0.5857864376
    };
    
    // Lower triangular coefficient matrix
    // C[i,j] = C[i*(i-1)/2 + j]
    static constexpr std::array<T, 1> C {
        -1.1715728753
    };

    // New step coefficients
    static constexpr std::array<T, S> M {
        0.8786796564, 
        0.2928932188,
    };

    // Error estimator coefficients
    static constexpr std::array<T, 2> E {
        0.2928932188, 
        0.2928932188
    };
    
    // Y_i ~= Y(T+H*Alpha[i])
    static constexpr std::array<T, 2> Alpha {
        0.0, 
        1.0
    };

    // Gamma_i = sum(Gamma[i,:])
    static constexpr std::array<T, 2> Gamma {
         1.7071067812, 
        -1.7071067812
    };
};


struct Ros4
{
    // Precision of integration parameters
    using T = double;

    // Estimator of local order
    static constexpr T ELO = 4.0;
    
    // Number of stages
    static constexpr int S = 4;

    // True: Evaluate f()
    // False: Reuse f() from previous stage
    static constexpr std::array<bool, S> EvalF { 
        true, 
        true, 
        true, 
        false 
    };

    // Lower triangular coefficient matrix
    // A[i,j] = A[i*(i-1)/2 + j]
    static constexpr std::array<T, 6> A { 
        2.0,
        1.867943637803922,
        0.2344449711399156,
        1.867943637803922,
        0.2344449711399156,
        0
    };
    
    // Lower triangular coefficient matrix
    // C[i,j] = C[i*(i-1)/2 + j]
    static constexpr std::array<T, 6> C { 
        -7.137615036412310,
         2.580708087951457,
         0.6515950076447975,
        -2.137148994382534,
        -0.3214669691237626,
        -0.6949742501781779
    };

    // New step coefficients
    static constexpr std::array<T, S> M { 
        2.255570073418735,
        0.2870493262186792,
        0.4353179431840180,
        1.093502252409163
    };

    // Error estimator coefficients
    static constexpr std::array<T, S> E { 
        -0.2815431932141155,
        -0.07276199124938920,
        -0.1082196201495311,
        -1.093502252409163
    };
    
    // Y_i ~= Y(T+H*Alpha[i])
    static constexpr std::array<T, S> Alpha { 
        0.0,
        1.145640000000000,
        0.6552168638155900,
        0.6552168638155900,
    };

    // Gamma_i = sum(Gamma[i,:])
    static constexpr std::array<T, S> Gamma { 
         0.5728200000000000,
        -1.769193891319233,
         0.7592633437920482,
        -0.1049021087100450
    };
};


template <size_t N,
          typename Parameters,
          typename LinearAlgebra>
class RosenbrockImpl
{
    using Vector = typename LinearAlgebra:: template Vector<N>;
    using Matrix = typename LinearAlgebra:: template Matrix<N, N>;
    using LUDecomp = typename LinearAlgebra:: template LUDecomp<Matrix>;

public:

    template <typename OdeFunction, typename OdeJacobian>
    static int integrate(
        OdeFunction fun,
        OdeJacobian jac,
        Vector & u,
        const double t0,
        const double tend,
        const double abstol=1.0,
        const double reltol=1e-3,
        const size_t maxstep=100000,
        const size_t maxdecomp=5,
        const double hmin=0,
        const double hmax=std::numeric_limits<double>::max(),
        const double hlim=1e-5,
        const double hfac=0.5,
        const double facmin=0.2,
        const double facmax=6,
        const double facrej=0.1,
        const double facsafe=0.9)
    {
        constexpr double eps = std::numeric_limits<double>::epsilon();
        constexpr double delmin = 1e-6;
        
        // Integration step count
        size_t nstep = 0;

        // Current integration time
        double t = t0;

        // Integration timestep
        double h = std::min(std::max(hmin, std::max(hmin, hlim)), std::max(hmax, tend-t0));

        // Stage vectors
        std::array<Vector, Parameters::S> K;

        // Time integration
        while (t <= tend) {
            if (nstep > maxstep) {
                // Error: too many steps
                return -1;
            }
            h = std::min(h, std::abs(tend-t));

            Vector f0;
            fun(f0, u, t);

            Matrix j0;
            jac(j0, u, t);

            // Finite difference approximation of df/dt
            const double tdel = std::sqrt(eps) * std::max(delmin, std::abs(t));
            const double tdel_inv = 1.0 / tdel;
            Vector dfdt;
            fun(dfdt, u, t+tdel);
            LinearAlgebra::aymx(dfdt, tdel_inv, f0);

            // Step calculation
            bool reject = false;
            while (true) {  
                // Construct step matrix: hgimj = 1/(h*gamma)*I - fJ(t)
                Matrix hgimj;
                LinearAlgebra::iama(hgimj, 1.0/(h*Parameters::Gamma[0]), j0);

                // Calculate LU decomposition of step matrix
                LUDecomp hgimj_decomp = LinearAlgebra::lu_decomposition(hgimj);
                
                // If decomposition is not invertable reduce step size and try again
                if (!LinearAlgebra::invertible(hgimj_decomp)) {         
                    double hbar = h * hfac;
                    for (size_t ndecomp=1; ndecomp<maxdecomp; ++ndecomp) {
                        LinearAlgebra::iama(hgimj, 1.0/(hbar*Parameters::Gamma[0]), j0);
                        LinearAlgebra::update_decomposition(hgimj_decomp, hgimj);
                        if (LinearAlgebra::invertible(hgimj_decomp)) {
                            // Decomp successful
                            break;
                        } else if (ndecomp == maxdecomp-1) {
                            // Error: Too many failed decompositions
                            return -1;
                        }
                        hbar *= hfac;
                    }
                    h = hbar;
                }

                // Calculate stages
                for (size_t stage=0; stage<Parameters::S; ++stage) {
                    Vector & sK = K[stage];
                    Vector fs;
                    if (stage == 0) {
                        LinearAlgebra::copy(fs, f0);
                    } else {
                        if (Parameters::EvalF[stage]) {
                            Vector ubar;
                            LinearAlgebra::copy(ubar, u);
                            for (size_t i=0; i<stage; ++i) {
                                double alpha = Parameters::A[stage*(stage-1)/2 + i];
                                LinearAlgebra::axpy(ubar, alpha, K[i]);
                            }
                            double tau = t + h*Parameters::Alpha[stage];
                            fun(fs, ubar, tau);
                        }
                    }
                    LinearAlgebra::copy(sK, fs);
                    for (size_t i=0; i<stage; ++i) {
                        double Ch = Parameters::C[stage*(stage-1)/2 + i] / h;
                        LinearAlgebra::axpy(sK, Ch, K[i]);
                    }
                    if (Parameters::Gamma[stage]) {
                        double hGamma = h * Parameters::Gamma[stage];
                        LinearAlgebra::axpy(sK, hGamma, dfdt);
                    }
                    LinearAlgebra::solve(hgimj_decomp, sK);
                } // Stage calculation

                // New solution
                Vector udot;
                LinearAlgebra::copy(udot, u);
                for (size_t i=0; i<Parameters::S; ++i) {
                    LinearAlgebra::axpy(udot, Parameters::M[i], K[i]);
                }

                // Estimate error
                Vector uerr;
                LinearAlgebra::scale(uerr, 0);
                for (size_t i=0; i<Parameters::S; ++i) {
                    LinearAlgebra::axpy(uerr, Parameters::E[i], K[i]);
                }
                double err = 0;
                for (size_t i=0; i<u.size(); ++i) {
                    double umax = std::max(std::abs(u[i]), std::abs(udot[i]));
                    double scale = abstol + reltol * umax;
                    err += (uerr[i] * uerr[i]) / (scale * scale);
                }
                err = std::sqrt(err / u.size());

                // New step size
                double hdot = h * std::min(facmax, std::max(facmin, facsafe/std::pow(err, 1.0/Parameters::ELO)));
                if ((err <= 1.0) || (h < hmin)) {
                    // Accept step and update solution
                    LinearAlgebra::copy(u, udot);
                    t += h;
                    // Set step size for next iteration
                    hdot = std::max(hmin, std::min(hdot, std::max(hmax, tend-t0)));
                    if (reject) {
                        h = std::min(hdot, h);
                    } else {
                        h = hdot;
                    }
                    // exit loop
                    break;
                }
                if (reject) {
                    // Consecutive rejections
                    // Back off step size by rejection factor and repeat step calculation
                    h *= facrej;
                } else {
                    // First rejection
                    // Retry step calculation with the new step size before backing off step size
                    reject = true;
                    h = hdot;
                }
            } // while (true)
        } // while (t <= tend)

        // Integration successful
        return 0;
    } // Integrate

}; // class Rosenbrock


} // namespace solver
