# really basic test suite just to make sure that the model part
# is working

using Test
using CytCKinetics

#=
    These tests are not written for correctness, but for whether the core code
    runs as expected. Correctness tests will be limited for those that rely
    only on reading information and likely not for manipulated information, though
    some exceptions may arise.
=#

#=
    These tests, unfortunately, run using the author's data and are not for release
    except under special exceptions (contact for more information).
=#

const configpath = joinpath(@__DIR__, "../data")
const configfile = joinpath(configpath, "WTC.toml")
const rawdatapath = joinpath(configpath, "WTC")
const setnames = ["Set1", "Set2"]

include("rawinput.jl")

@test 1+1 == 2
