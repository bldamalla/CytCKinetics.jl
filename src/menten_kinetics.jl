# By Bon Leif AMALLA (B4, 2022)

# Written for basic activity measurements. In using this initial version, there will
# be certain restrictions on the format of the input files and their filenames.
# More details can be found in the README file.

# Utility functions for fitting using LsqFit

using StaticArrays
export FitResults

"""
    FitResults{T<:AbstractFloat}

Convinience container for storing results from fitting best fit line to single
series data from a Series set.
"""
struct FitResults{T<:AbstractFloat}
    slope::T
    N::Int
    r2::T
end

"""
    menten(x, params)
    menten_inplace(F, x, params)

Michaelis-Menten model function under constant enzyme concentration. The initial
substrate concentration is `x`. `params` is an iterable with two elements, where the
first element is ``V_{max}`` and the second element is ``K_M``.

```math
v_0(x) = \\dfrac{V_{max} \\cdot x}{K_{M} + x}
```

Note: Usually `params` is a `Vector{T<:Real}` is passed to `LsqFit.curve_fit` 
during fitting.  Inplace function is provided for fitting when needed.

See also: [`menten_jac`](@ref)
"""
menten(s, p) = @. (p[1] * s) / (p[2] + s)
@doc (@doc menten) menten_inplace(F, s, p) = (@. F = (p[1] * s) / (p[2] + s))

"""
    menten_jac(x, params)
    menten_jac_inplace(J, x, params)

Jacobian of the Michaelis-Menten model function under constant enzyme concentration. The
initial substrate concentration is `x`. `params` is an iterable with two elements, where
the first element is ``V_{max}`` and the second element is ``K_M``.

Note: Usually `params` is a `Vector{T<:Real}` is passed to `LsqFit.curve_fit` 
during fitting.  Inplace function is provided for fitting when needed.

See also: [`menten`](@ref)
"""
function menten_jac(s, p)
    J = Matrix{Float64}(undef, length(s), length(p))
    @. J[:,1] = s / (s + p[2])
    @. @views J[:,2] = -p[1] * J[:,1] / (s + p[2])
    J
end

## this is mostly for dispatching by the ribbon functions when plotting
function menten_jac(s::Number, p)
    ## do the jacobian stuff here
    Ja = s / (s + p[2])
    Jb = -p[1] * Ja / (s + p[2])
    return @SMatrix [Ja Jb]
end

function menten_jac_inplace(J::Array{Float64, 2}, s, p)
    @. J[:,1] = s / (s + p[2])
    @. @views J[:,2] = -p[1] * J[:,1] / (s + p[2])
end
@doc (@doc menten_jac) menten_jac_inplace
