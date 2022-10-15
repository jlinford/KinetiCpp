#include <cmath>
#include <array>

struct Ros2
{
    // Precision of integration parameters
    using T = double;

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
