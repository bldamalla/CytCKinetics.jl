# By Bon Leif AMALLA (B4, 2022)

# Written for basic activity measurements. In using this initial version, there will
# be certain restrictions on the format of the input files and their filenames.
# More details can be found in the README file.

## SERIES SET TYPES AND FITTING

# probably a no op anyway --- just to make sure
using StaticArrays
using LsqFit

export AbsorbanceSeries, SeriesSet, SeriesSetResults
export blankedframe, rateconst
export fit

"""
    AbsorbanceSeries

Container for information regarding how to control absorbance data
for kinetics manipulations. Specifically developed for cytochrome c
kinetics.
"""
struct AbsorbanceSeries{T<:AbstractFloat}
    ardata::AbsorbanceRaw

    # information needed for processing
    framestart::T
    framestop::T
    blanktime::T
end

"""
    SeriesSet{T<:Real}

Representation of a set of `AbsorbanceSeries` objects and other
corresponding information.

Note: `menten=true` means series data has to be fit using nonlinear
least squares regression.
"""
struct SeriesSet{T<:Real}
    series::Vector{AbsorbanceSeries}
    concentrations::Vector{T}
    menten::Bool
    
    # other metadata contained in the TOML file
    metadata::Dict{String,Any}
end

"""
    SeriesSetResults{T<:Real}

Results obtained after calculations of initial reaction rates and
Michaelis-Menten parameters from nonlinear least squares regression.

Note: `SeriesSets` with `menten=false` have `fitparams`, `stderrors`,
and `covmatrix` fields set to `nothing`.
"""
struct SeriesSetResults{T<:Real}
    # plotting essentials
    concentrations::Vector{T}
    initrates::Vector{T}

    # results from LsqFit regression
    fitparams::Union{Nothing,SVector{2,T}}
    stderrors::Union{Nothing,Tuple{T,T}}
    covmatrix::Union{Nothing,SMatrix{2,2,T}}
end

"""
    blankedframe(ser::AbsorbanceSeries)

Return an `ARManip` containing blanked absorbance data from series
information.
"""
function blankedframe(ser::AbsorbanceSeries)
    (; ardata, framestart, framestop) = ser
    arview = getARView(ardata, framestart, framestop)

    (; blanktime) = ser
    blankidx = time2index(ardata, blanktime)
    
    return arview - ardata[blankidx]
end

"""
    rateconst(ser::AbsorbanceSeries; r2thresh, minthresh)

Simplified routine to calculate initial rate constant for cytochrome
_c_ oxidation, almost specifically (purpose of the module, really).
This routine should also work for any first order reactions, or
basically any routine involving an approximate exponential decrease
in absorbance value.

Note: Also returns the number of points used in fitting and the
coefficient of deterination ``r^2`` from fitting.
"""
function rateconst(ser::AbsorbanceSeries; r2thresh=0.96, minthresh=20)
    logAbs = blankedframe(ser) |> log

    # find appropriate N and r2 for linear fitting
    # based on needed linearity parameters
    N, r2fitted = thresholdfit(logAbs; r2thresh=r2thresh, minthresh=minthresh)

    # actual linear fitting
    k, _ = linfit(logAbs.times[1:N], logAbs.values[1:N])

    return FitResults(k, N, r2fitted)
end

"""
    fit(SeriesSetResults, serset::SeriesSet, fitstart; r2thresh, minthresh)

Obtain initial reaction rates and concentrations for plotting. If the
series set is not admissible for fitting, _i.e._ `serset.menten=false`
then no fitting will be done. Instead, just rates and concentrations
will be returned.

Note: `fitstart` is an optional positional argument. This specifies the starting vector
for fitting with the Michaelis-Menten model.
"""
function fit(::Type{SeriesSetResults}, serset::SeriesSet,
             fitstart=[60,0.4]; kwargs...)
    (; series, concentrations) = serset
    initrates = map(zip(series, concentrations)) do (ser, conc)
        k = rateconst(ser; kwargs...)
        conc * -k.slope * 1000
    end

    (; menten) = serset
    menten || return SeriesSetResults(
        concentrations, initrates,
        # then no fitting is done
        nothing, nothing, nothing
    )

    # start the fitting routine
    fitmodel = curve_fit(menten_inplace, menten_jac_inplace,
                         concentrations, initrates, fitstart;
                         inplace=true)

    fitmodel.converged || @warn "Fitting did not converge..."

    fitparams = fitmodel.param
    stderrors = stderror(fitmodel)
    covmat = estimate_covar(fitmodel)
    return SeriesSetResults(
        concentrations, initrates,
        fitparams, stderrors, covmat
    )
end
