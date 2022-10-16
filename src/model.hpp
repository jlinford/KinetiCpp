#pragma once
#include <array>
#include "rosenbrock.hpp"



template <typename T, typename Method, typename Mechanism>
struct Model
{
    using var_t = std::array<T, Mechanism::NVAR>;
    using fix_t = std::array<T, Mechanism::NFIX>;

    size_t const maxstep = 100000;
    double const facmin = 0.2;
    double const facmax = 6;
    double const facrej = 0.1;
    double const facsafe = 0.9;

    var_t ode_function(double t, var_t const & var, fix_t const & fix)
    {
        return var;
#if 0
    // Compute equation rates
    for(j=0; j<EqnNr; j++) {
        used = 0;
        for (i = 0; i < VarNr; i++) {
            if ( Stoich[i][j] != 0 ) { 
                used = 1;
                break;
            }
        }
        if ( used ) {    
            prod = RConst( j );
            for (i = 0; i < VarNr; i++) {
                for (k = 1; k <= (int)Stoich_Left[i][j]; k++ ) {
                    prod = Mul( prod, Elm( V, i ) ); 
                }
            }
            for ( ; i < SpcNr; i++) {
                for (k = 1; k <= (int)Stoich_Left[i][j]; k++ ) {
                    prod = Mul( prod, Elm( F, i - VarNr ) );
                }
            }
            Assign( Elm( A, j ), prod );
        }
    }

    // Aggregate function
    for (i = 0; i < VarNr; i++) {
        sum = Const(0);
        for (j = 0; j < EqnNr; j++) 
        sum = Add( sum, Mul( Const( Stoich[i][j] ), Elm( A, j ) ) );
        Assign( Elm( Vdot, i ), sum );
    }
#endif
    }

    int integrate(
        double t0,
        double tend,
        double tdel,
        var_t & var,
        fix_t const & fix,
        double abstol=1.0,
        double reltol=1e-3,
        double hmin=0,
        double hmax=std::numeric_limits<double>::max(),
        double hstart=0,
        double hlim=1e-5)
    {
        using std::min;
        using std::max;
        using std::abs;

        double constexpr eps = std::numeric_limits<double>::epsilon();

        // Integration step count
        size_t nstep = 0;

        // Current integration time
        double t = t0;

        // Integration timestep
        double h;
        hmax = max(hmax, tend-t0);
        hstart = max(hmin, hlim);
        h = min(max(hmin, hstart), hmax);
        if (h <= 10*tdel) {
            h = hlim;
        }

        bool reject_last = false;
        bool reject_more = false;

        while (t <= tend) {
            if (nstep > maxstep) {
                // Error: too many steps
                return -1;
            }
            h = min(h, abs(tend-t));

            auto f0 = ode_function(t, var, fix);
        //     if (time_dependent) {
        //         dFdT = fdt(t, tdel, y, f0);
        //     }

        //     jac0 = jac(t, y);

        //     while (true) {
        //         //CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), Jac0,Ghimj,Pivot,Singular)
        //         if (singular) {
        //             // Error: too many failed decompositions
        //             return -1;
        //         }

        //         for (int istage=0; i<Method::S; ++i) {
        //             int offset = istage * N;
        //             if (istage == 0) {
        //                 // Fcn(1:N) = Fcn0(1:N)
        //             } else if Method::EvalF(istage) {
        //                 // Ynew(1:N) = Y(1:N)
        //                 for (int j=0; j<istage; ++j) {
        //                     //CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), K(N*(j-1)+1),1,Ynew,1)
        //                 }
        //                 // Tau = T + ros_Alpha(istage)*Direction*H
        //                 // CALL FunTemplate(Tau,Ynew,Fcn)
        //             }

        //             // K(ioffset+1:ioffset+N) = Fcn(1:N)
        //             for (int j=0; j<istage; ++j) {
        //                 //HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
        //                 //CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
        //             }
        //             if (time_dependent && Method:Gamma[istage] != 0) {
        //                 // HG = Direction*H*ros_Gamma(istage)
        //                 //CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
        //             }
        //             // CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))
        //         } // for(i...)

        //         ynew = y;
        //         for (j=0; j<Method::S; ++j) {
        //             // CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
        //         }

        //         yerr = 0;
        //         for (j=0; j<Method::S; ++j) {
        //             // CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
        //         }
        //         // Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

        //         // New step size is bounded by FacMin <= Hnew/H <= FacMax
        //         //Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
        //         hnew = h*fac;

        //         // Adjust step size according to error magnitude
        //         if ((err <= 1.0) || (h <= hmin)) {
        //             // Accept step
        //             y = ynew;
        //             t += h;
        //             hnew = max(hmin,min(hnew,hmax));
        //             if (reject_last) {
        //                 hnew = min(hnew, h);
        //             }
        //             RSTATUS(Nhexit) = H
        //             RSTATUS(Nhnew)  = Hnew
        //             RSTATUS(Ntexit) = T
        //             RejectLastH = .FALSE.
        //             RejectMoreH = .FALSE.
        //             H = Hnew
        //             // Exit while(true) loop
        //             break;
        //         } else {
        //             // Reject step
        //             if (reject_more) {
        //                 hnew = h*facrej;
        //             }
        //             reject_more = reject_last;
        //             reject_last = true;
        //             h = hnew;
        //         }
        //     } // while (true)
        } // while (t <= tend)

        // Integration successful
        return 0;
    }
};