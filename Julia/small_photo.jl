#
# A translation of the "Small" mechanism from KPP (Sandu et al.)
#
# <R1>  O2   + hv = 2O       : (2.643E-10) * SUN*SUN*SUN;
# <R2>  O    + O2 = O3       : (8.018E-17);
# <R3>  O3   + hv = O   + O2 : (6.120E-04) * SUN;
# <R4>  O    + O3 = 2O2      : (1.576E-15);
# <R5>  O3   + hv = O1D + O2 : (1.070E-03) * SUN*SUN;
# <R6>  O1D  + M  = O   + M  : (7.110E-11);
# <R7>  O1D  + O3 = 2O2      : (1.200E-10);
# <R8>  NO   + O3 = NO2 + O2 : (6.062E-15);
# <R9>  NO2  + O  = NO  + O2 : (1.069E-11);
# <R10> NO2  + hv = NO  + O  : (1.289E-02) * SUN;
# 
# John Linford <jlinford@redhpc.com>
#
# Copyright 2022 John Linford
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
#   1. Redistributions of source code must retain the above copyright 
#      notice, this list of conditions and the following disclaimer.
#   2. Redistributions in binary form must reproduce the above copyright 
#      notice, this list of conditions and the following disclaimer in the 
#      documentation and/or other materials provided with the distribution.
#   3. Neither the name of the copyright holder nor the names of its 
#      contributors may be used to endorse or promote products derived from
#      this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.
#

using LinearAlgebra
using DifferentialEquations
using SparseArrays
using Plots

# Reaction equations and rates.
# Equations are declared as tuples of `(Reactants, Products, Rate)`.
# - Reactants and products are declared as lists of equation terms.  
#   The terms will summed in parsing.
# - Equation terms are declared as `Symbol` species names as with an 
#   optional `Float64` term coefficient.  If omitted, the coefficient 
#   is assumed to be 1.
# - Reaction rates may be a `Float64` value, or a callable `rate = f(x)`.
reactions = [
    ([:O2],       [2, :O],     sun -> (2.643e-10 * sun^3)),
    ([:O, :O2],   [:O3],       8.018e-17),
    ([:O3],       [:O, :O2],   sun -> (6.12e-04 * sun)),
    ([:O, :O3],   [2, :O2],    1.576e-15),
    ([:O3],       [:O1D, :O2], sun -> (1.07e-03 * sun^2)),
    ([:O1D, :M],  [:O, :M],    7.11e-11),
    ([:O1D, :O3], [2, :O2],    1.2e-10),
    ([:NO, :O3],  [:NO2, :O2], 6.062e-15),
    ([:NO2, :O],  [:NO, :O2],  1.069e-11),
    ([:NO2],      [:NO, :O],   sun -> (1.289e-02 * sun))
]

# Variable species concentrations change according to the law of mass action kinetics
var_spec = [
    :O1D, 
    :O, 
    :O3, 
    :NO, 
    :NO2
]

# Fixed species concentrations are determined by physical factors
fix_spec = [
    :M, 
    :O2
]

# Bind matrix column numbers to species symbols
# The order of `spec` determines the sparsity pattern of the stoichiometric matrix
# Fixed species must follow variable species in `spec`
spec = [var_spec ; fix_spec]
for (i, s) in enumerate(spec)
    eval(:($s = $i))
end

# Convenience variables
nvar = length(var_spec)
nfix = length(fix_spec)
nspec = length(spec)
nreact = length(reactions)

# Stoichiometric matrics
# LHS and RHS stoichiometric matrices are stored transposed to improve data locality
# Aggregate stoichiometric matrix only stores variable species
lhs_stoich = spzeros((nspec, nreact))
rhs_stoich = spzeros((nspec, nreact))
agg_stoich = spzeros((nreact, nvar))

# Parse the reactions to build the stoichiometric matrices
let
    function parse_terms(terms)
        coef, spec = nothing, nothing
        Channel() do channel
            for t in terms
                if t isa Number
                    if !(isnothing(coef) && isnothing(spec))
                        error("Invalid reaction: unexpected number: ", terms)
                    end
                    coef = Float64(t)
                    spec = nothing
                    continue
                elseif t isa Symbol
                    if !isnothing(spec)
                        error("Invalid reaction: unexpected string: ", terms)
                    end
                    spec = t
                else
                    error("Invalid reaction: invalid term: ", terms)
                end
                if isnothing(coef)
                    coef = 1.0
                end
                put!(channel, (coef, spec))
                coef, spec = nothing, nothing
            end
        end 
    end
    for (i, rct) in enumerate(reactions)
        lhs, rhs, _ = rct
        for (coef, spec) in parse_terms(lhs)
            j = eval(spec)
            lhs_stoich[j,i] += coef
            if j <= nvar
                agg_stoich[i,j] -= coef
            end
        end
        for (coef, spec) in parse_terms(rhs)
            j = eval(spec)
            rhs_stoich[j,i] += coef
            if j <= nvar
                agg_stoich[i,j] += coef
            end
        end
    end
end 

# Calculate Jacobian sparsity structure
structJ = Matrix{Bool}(I, nspec, nspec)
let
    for i in 1:nvar, j in 1:nspec, k in 1:nreact
        if agg_stoich[k,i]*lhs_stoich[j,k] != 0
            structJ[i,j] = true
        end
    end
end


"""
    sunlight(t)

Returns sunlight intensity between 0 and 1.0.

# Arguments 
- `t`: Time in seconds
"""
function sunlight(t)
    sunrise = 4.5 * 3600    # 4:30am
    sunset  = 19.5 * 3600   # 7:30pm
    if (t < sunrise) || (t > sunset)
        return 0
    end
    θ = abs((2*t-sunrise-sunset)/(sunset-sunrise))^2
    return (1 + cos(π*θ)) / 2
end
# Sanity check
#plot([sunlight(t) for t::Float64 in 0:60:24*3600])

"""
    calc_rates(T, t)

Returns an array of reaction rate constants at time `t`.
# Arguments
- `T`: the type of the rate constant.  Necessary to allow Dual types to propogate from DifferentialEquations.jl
- `t`: Time in seconds
"""
function calc_rates(T, t)::AbstractArray{T}
    sun = sunlight(t)
    rates = zeros(T, nreact)
    for (i, rct) in enumerate(reactions)
        _, _, rate = rct
        if rate isa Number
            rates[i] = rate
        else
            rates[i] = rate(sun)
        end
    end
    return rates
end
# Sanity check
#plot([calc_rates(Float64, t)[10] for t::Float64 in 0:60:24*3600])

"""
    f!(du, u, p, t)

Aggregate production/destruction function.

`T` is captured to allow Dual types to propogate from DifferentialEquations.jl
"""
function f!(du::AbstractVector{T}, u, p, t) where T
    rates = calc_rates(T, t)
    rate_prod = similar(rates)
    for i in eachindex(rates)
        rate_prod[i] = rates[i] * prod(u.^lhs_stoich[:,i])
    end
    for j in 1:nvar
        du[j] = sum(agg_stoich[:,j] .* rate_prod[:])
    end
    # Fixed species concentrations do not change
    du[nvar+1:nspec] .= 0
end

"""
    fJ!(J, u, p, t)

Jacobian of the ODE function.

`T` is captured to allow Dual types to propogate from DifferentialEquations.jl
"""
function fJ!(J::AbstractArray{T}, u, p, t) where T
    rates = calc_rates(T, t)
    B = spzeros(T, (nreact, nvar))
    for i in 1:nreact, j in 1:nvar
        if lhs_stoich[j,i] != 0
            p = rates[i]
            p *= prod(u[1:j-1].^lhs_stoich[1:j-1,i])
            p *= lhs_stoich[j,i] * u[j]^(lhs_stoich[j,i]-1)
            p *= prod(u[j+1:nspec].^lhs_stoich[j+1:nspec,i])
            B[i,j] = p
        end
    end
    for i in 1:nvar, j in 1:nvar
        if structJ[i,j]
            J[i,j] = sum(agg_stoich[:,i] .* B[:,j])
        end
    end
end


# Time integration
let
    # Initial concentrations
    # Order must match `spec`
    u0 = [
        # Variable species concentrations
        9.906e+01,
        6.624e+08,
        5.326e+11,
        8.725e+08,
        2.240e+08,
        # Fixed species concentrations
        8.120e+16,
        1.697e+16
    ]

    # Integration time in seconds
    t0 = 0 * 3600
    tend = 24 * 3600

    # Integration tolerances
    abstol = 1.0
    reltol = 1e-3

    # Define the ODE problem.
    # Map species names to elements of the solution vector
    ff = ODEFunction(f! ; jac=fJ!, jac_prototype=float.(structJ), syms=spec)

    # Solve the ODE problem
    sol = solve(ODEProblem(ff, u0, (t0, tend)), 
                Rosenbrock23(), 
                abstol=abstol, 
                reltol=reltol)

    # Plot the solution
    plot(plot(sol, idxs=[2, 4, 5]),
         plot(sol, idxs=[3]),
         plot(sol, idxs=[1, 6, 7]))
end
