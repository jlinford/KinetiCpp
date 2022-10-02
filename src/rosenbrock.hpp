#include <cmath>
#include <array>

using std::size_t;

template <typename T>
struct Ros2
{
    // Partial gamma
    static constexpr T G = 1.0 + 1.0/std::sqrt(2.0);

    // Estimator of local order
    static constexpr T ELO = 2.0;
    
    // Number of stages
    static constexpr int S = 2;

    // True: Evaluate f()
    // False: Reuse f() from previous stage
    static constexpr std::array<bool, S> EvalF {true, true};

    // Lower triangular coefficient matrix
    // A[i,j] = A[i*(i-1)/2 + j]
    static constexpr std::array<T, 1> A {1.0/G};
    
    // Lower triangular coefficient matrix
    // C[i,j] = C[i*(i-1)/2 + j]
    static constexpr std::array<T, 1> C {-2.0/G};

    // New step coefficients
    static constexpr std::array<T, 2> M {3.0/(2.0*G), 1.0/(2.0*G)};

    // Error estimator coefficients
    static constexpr std::array<T, 2> E {1.0/(2.0*G), 1.0/(2.0*G)};
    
    // Y_i ~= Y(T+H*Alpha[i])
    static constexpr std::array<T, 2> Alpha {0, 1};

    // Gamma_i = sum(Gamma[i,:])
    static constexpr std::array<T, 2> Gamma {G, -G};
};


template <typename Method, typename T, size_t N>
struct Rosenbrock
{
    using elem_t = T;
    using vector_t = std::array<T, N>;

    Rosenbrock(bool autonomous=true,
               size_t maxstep=10000,
               double hmin=0,
               double hmax=0,
               double hstart=0,
               double facmin=0.2,
               double facmax=6,
               double facrej=0.1,
               double facsafe=0.9)
    {

    }

    int Integrate(
        double && t0, double && tend, 
        vector_t & y, 
        vector_t const & abstol,
        vector_t const & reltol,
        bool time_dependent=true)
    {
        using std::min;
        using std::max;
        using std::abs;

        constexpr double delmin = 1.0e-5;

        double t = t0;
        
        double h = min(max(abs(hmin),abs(hstart)), abs(hmax));
        if (abs(h) <= 10*tdel) {
            h = delmin;
        }

        bool reject_last = false;
        bool reject_more = false;

        while (t <= tend) {
            if (nstep > maxstep) {
                // Error: too many steps
                return -1;
            }
            h = min(h, abs(tend-t));

            f(t, y, f0);
            if (time_dependent) {
                dFdT = fdt(t, tdel, y, f0);
            }

            jac0 = jac(t, y);

            while (true) {
                //CALL ros_PrepareMatrix(H,Direction,ros_Gamma(1), Jac0,Ghimj,Pivot,Singular)
                if (singular) {
                    // Error: too many failed decompositions
                    return -1;
                }

                for (int istage=0; i<Method::S; ++i) {
                    int offset = istage * N;
                    if (istage == 0) {
                        // Fcn(1:N) = Fcn0(1:N)
                    } else if Method::EvalF(istage) {
                        // Ynew(1:N) = Y(1:N)
                        for (int j=0; j<istage; ++j) {
                            //CALL WAXPY(N,ros_A((istage-1)*(istage-2)/2+j), K(N*(j-1)+1),1,Ynew,1)
                        }
                        // Tau = T + ros_Alpha(istage)*Direction*H
                        // CALL FunTemplate(Tau,Ynew,Fcn)
                    }

                    // K(ioffset+1:ioffset+N) = Fcn(1:N)
                    for (int j=0; j<istage; ++j) {
                        //HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
                        //CALL WAXPY(N,HC,K(N*(j-1)+1),1,K(ioffset+1),1)
                    }
                    if (time_dependent && Method:Gamma[istage] != 0) {
                        // HG = Direction*H*ros_Gamma(istage)
                        //CALL WAXPY(N,HG,dFdT,1,K(ioffset+1),1)
                    }
                    // CALL ros_Solve(Ghimj, Pivot, K(ioffset+1))
                } // for(i...)

                ynew = y;
                for (j=0; j<Method::S; ++j) {
                    // CALL WAXPY(N,ros_M(j),K(N*(j-1)+1),1,Ynew,1)
                }

                yerr = 0;
                for (j=0; j<Method::S; ++j) {
                    // CALL WAXPY(N,ros_E(j),K(N*(j-1)+1),1,Yerr,1)
                }
                // Err = ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )

                // New step size is bounded by FacMin <= Hnew/H <= FacMax
                //Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
                hnew = h*fac;

                // Adjust step size according to error magnitude
                if ((err <= 1.0) || (h <= hmin)) {
                    // Accept step
                    y = ynew;
                    t += h;
                    hnew = max(hmin,min(hnew,hmax));
                    if (reject_last) {
                        hnew = min(hnew, h);
                    }
                    RSTATUS(Nhexit) = H
                    RSTATUS(Nhnew)  = Hnew
                    RSTATUS(Ntexit) = T
                    RejectLastH = .FALSE.
                    RejectMoreH = .FALSE.
                    H = Hnew
                    // Exit while(true) loop
                    break;
                } else {
                    // Reject step
                    if (reject_more) {
                        hnew = h*facrej;
                    }
                    reject_more = reject_last;
                    reject_last = true;
                    h = hnew;
                }
            } // while (true)
        } // while (t <= tend)
        
        // Integration successful
        return 0;
    }

    size_t maxstep = 10000;
    double hmin = 0;
    double hmax = 0;
    double hstart = 0;
    double facmin = 0.2;
    double facmax = 6;
    double facrej = 0.1;
    double facsafe = 0.9;

};