# By Bon Leif AMALLA (B4, 2022)

# Written for basic activity measurements. In using this initial version, there will
# be certain restrictions on the format of the input files and their filenames.
# More details can be found in the README file.

## MERGE SERIES SETS INTO SET GROUPS FOR CLEANER STATISTICAL TREATMENT
export SetGroup, SetGroupResults

struct SetGroup{T<:Real}
    sets::Vector{SeriesSet}
    # menten::Bool should be implicitly true otherwise there'd be no other
    # reason to use groups

    # other metadata contained in the TOML file
    metadata::Dict{String,Any}
end

struct SetGroupResults{T<:Real}
    # plotting essentials
    concentrations::Vector{T}
    initrates::Vector{T}

    # results from LsqFit regression
    fitparams::SVector{2,T}
    stderrors::Tuple{T,T}
    covmatrix::SMatrix{2,2,T}
end

function setmerge(::Type{SetGroup}, serset::SeriesSet{T}...) where T
    sets = collect(serset)
    emptydict = Dict{String,Any}()
    return SetGroup{T}(sets, emptydict)
end

"""
    fit(SetGroupResults, setgrp::SetGroup, fitstart; r2thresh, minthresh)

Fit combined series set initial reaction rates and concentrations using the
Michaelis-Menten model. In contrast to the `fit` method for `SeriesSetResults`, this
method assumes it to be true.

Note: `fitstart` is an optional position argument. This specifies the starting vector
for fitting with the Michaelis-Menten model. Internally uses the inplace functions.
"""
function fit(::Type{SetGroupResults}, setgrp::SetGroup,
             fitstart=[60,0.4]; kwargs...)
    # appropriately merge the concentrations...
    concentrations = reduce(vcat, (ser.concentrations for ser in setgrp.sets))
    # ... and initial reaction rates
    initrates = reduce(vcat, (initialrates(ser; kwargs...) for ser in setgrp.sets))

    fitmodel = curve_fit(menten_inplace, menten_jac_inplace,
                         concentrations, initrates, fitstart;
                         inplace=true)

    fitmodel.converged || @warn "Fitting did not converge..."

    fitparams = fitmodel.param |> SVector{2}
    stderrors = stderror(fitmodel) |> q->tuple(q...)
    covmat = estimate_covar(fitmodel) |> SMatrix{2,2}

    return SetGroupResults(
        concentrations, initrates,
        fitparams, stderrors, covmat
    )
end
