## commons file for scripts

import TOML

const datapath = joinpath(@__DIR__, "../data/")
const resultspath = joinpath(datapath, "results.toml")
parse_setargument(str) = split(str, ':')

"""
    SeriesData{T<:Real}

Convenience struct for containing the information needed to plot
at the end of analyses.
"""
struct SeriesData{T<:Real}
    ## these are the actual data points
    concentrations::Vector{T}
    initialrates::Vector{T}

    # fitting from using LsqFit
    mentenparams::Union{Tuple{T,T},Nothing}
    ## uncomment this when there are plans to make confidence bands for multiset or
    ## multigroup routines
    # covmatrix::Union{Matrix{T},Nothing}

    # these are the properties of the data file studied
    # dates of measurement likely arent needed for plotting for
    # the meantime, so it may be okay to omit those for now
    strain::String
    setname::String
end

function getmeta(strain)
    open(joinpath(datapath, "./$(strain).toml")) do io
        TOML.parse(io)
    end
end

function getsetroot(strainmeta, set)
    joinpath(datapath, strainmeta["rootfolder"], strainmeta[set]["subfolder"])
end

getsetseries(strainmeta, set) = strainmeta[set]["Series"]

function threadmap(f::Function, iter)
    tasks = map(iter) do item
        Threads.@spawn f(item)
    end
    fetch.(tasks)
end

using StaticArrays
using CytCKinetics

const fitstart = [60, 0.5]

@info "Loaded commons scripts..."
