# KinetiCpp

KinetiCpp is a system for simulating mass action kinetics implemented as pure C++23 headers.  

 * Compile-time formulation of stochiometric matrices.  
 * Extremely efficient implementations of the ODE function and Jacobian.
 * A Rosenbrock time-stepping integrator well-suited for solving highly stiff ODE systems.
 * A user-friendly domain specific language for defining chemical reaction networks.

Here's a functional example of a Chapman-like mechanism implemented with KinetiCpp:

```c++
using Chapman = Mechanism <
    VariableSpecies {
        O1D || O, 
        O1  || O, 
        O3  || O*3, 
        NO  || N + O,
        NO2 || N + O*2
    },
    FixedSpecies {
        M  || N*2 + O*2, 
        O2 || O*2
    },
    O2       >= 2 * O1   || [](double t) { return 2.643e-10 * std::pow(sunlight(t), 3); },
    O1 + O2  >= O3       || 8.018e-17,
    O3       >= O1 + O2  || [](double t) { return 6.12e-04 * sunlight(t); },
    O1 + O3  >= 2 * O2   || 1.576e-15,
    O3       >= O1D + O2 || [](double t) { return 1.07e-03 * std::pow(sunlight(t), 2); },
    O1D + M  >= O1 + M   || 7.11e-11, 
    O1D + O3 >= 2 * O2   || 1.2e-10,
    NO + O3  >= NO2 + O2 || 6.062e-15, 
    NO2 + O1 >= NO + O2  || 1.069e-11,
    NO2      >= NO + O1  || [](double t) { return 1.289e-02 * sunlight(t); }
>;
```

See [examples](examples) for more details.

----

John Linford <jlinford@redhpc.com>

