# setmanips/basics.jl -- handling set information

# probably a no op anyway --- just to make sure
using StaticArrays

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

struct SeriesSet{T<:Real}
    series::Vector{AbsorbanceSeries}
    concentrations::Vector{T}
    menten::Bool
    
    # other metadata contained in the TOML file
    metadata::Dict{String,Any}
end

struct SeriesSetResults{T<:Real}
    # plotting essentials
    concentrations::Vector{T}
    initrates::Vector{T}

    # results from LsqFit regression
    fitparams::Union{Nothing,SVector{2,T}}
    stderrors::Union{Nothing,Tuple{T,T}}
    covmatrix::Union{Nothing,SMatrix{2,2,T}}
end
