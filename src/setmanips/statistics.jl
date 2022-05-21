# By Bon Leif AMALLA (B4, 2022)

# Written for basic activity measurements. In using this initial version, there will
# be certain restrictions on the format of the input files and their filenames.
# More details can be found in the README file.

## MORE FORMAL TREATMENT OF SERIES SET AND GROUP SET STATISTICS

"""
    fitstderror(ssr, x::Number)
    fse(ssr, x)

Represents the standard error of the dependent variable, _i.e._ initial reaction
rates, at the given initial substrate concentration `x`. Calculated as follows:

```math
    \\sigma_y = J'(x, p_0) \\cdot G(p_0) \\cdot J(x, p_0),
```

where ``J(x, p_0)`` is the Jacobian of the Michaelis-Menten model equation evaluated
at the given substrate concentration and the fitted parameters from regression. Its
transpose is the primed term. ``G`(p_0)`` is the covariance matrix calculated by
`LsqFit`.

Note: `fse` is an alias to `fitstderror`.
"""
function fitstderror(ssr, s::Number)
    ∇ = menten_jac(s, ssr.fitparams)
    variance = ∇ * ssr.covmatrix * transpose(∇)
    return sqrt(variance)
end
const fse = fitstderror
@doc (@doc fitstderror) fse