# By Bon Leif AMALLA (B4, 2022)

# Written for basic activity measurements. In using this initial version, there will
# be certain restrictions on the format of the input files and their filenames.
# More details can be found in the README file.

## TYPE DEFINITIONS AND MANIPULATION -- INCLUDES READING EXPECTED TXT FORMAT

export AbsorbanceRaw, ARView
export starttime, stoptime, timestep
export time2index
export getARView

const FT = float(Int)

"""
    AbsorbanceRaw{T} <: AbstractVector{Tuple{T,T}}

Struct containing raw absorbance data obtained from time series steady-state
kinetic experiments. As of writing, this is the result of reading output files
from the used spectrometer.
"""
struct AbsorbanceRaw{T} <: AbstractVector{Tuple{T,T}}
    times::Vector{T}
    values::Vector{T}
end
const AR = AbsorbanceRaw
function AbsorbanceRaw(times, values)
    @assert length(times) == length(values)
    T = promote_type(eltype(times), eltype(values))
    return AbsorbanceRaw{T}(times, values)
end
Base.size(ar::AR) = (length(ar.times),)

function Base.getindex(ar::AR, i::Int)
    @boundscheck checkbounds(ar, i)
    return @inbounds (ar.times[i], ar.values[i])
end

starttime(ar::AR) = first(ar.times)
stoptime(ar::AR) = last(ar.times)
timestep(ar::AR) = @inbounds ar.times[begin+1] - first(ar.times)

time2index(ar::AR, t::Real) = floor(Int, (t - starttime(ar)) / timestep(ar)) + 1

"""
    ARView{T} <: AbstractVector{T}

Lazy reference to an `AbsorbanceRaw` object. Can also be referenced similar to a
vector. Note that the indices of the `ARView` are still the same used by the
parent `AbsorbanceRaw` object.
"""
struct ARView{T} <: AbstractVector{Tuple{T,T}}
    ar::AR{T}
    range::UnitRange{Int}

    function ARView(ar::AR{T}, start, stop) where T
        checkbounds(ar, start)
        checkbounds(ar, stop)
        @assert stop > start "start index must be strictly less than stop index"
        return new{T}(ar, start:stop)
    end
end
Base.size(arv::ARView) = size(arv.range)
Base.axes(arv::ARView) = (arv.range,)

Base.checkbounds(::Type{Bool}, arv::ARView, i::Int) = i ∈ arv.range
Base.checkbounds(arv::ARView, i::Int) = begin
    checkbounds(Bool, arv, i) || throw(BoundsError(arv, i))
    nothing
end

function Base.getindex(arv::ARView, i::Int)
    @boundscheck checkbounds(arv, i)
    return @inbounds arv.ar[i]
end

starttime(arv::ARView) = @inbounds arv.ar.times[first(arv.range)]
stoptime(arv::ARView) = @inbounds arv.ar.times[last(arv.range)]
timestep(arv::ARView) = timestep(arv.ar)

time2index(arv::ARView, t::Real) = time2index(arv.ar, t)

function Base.getproperty(arv::ARView, s::Symbol)
    if s == :times
        return arv.ar.times[arv.range]
    elseif s == :values
        return arv.ar.values[arv.range]
    else
        Base.getfield(arv, s)
    end
end

function Base.read(io::IO, ::Type{AbsorbanceRaw}, assert_=false)
    ## do stuff here
    parsestart = false
    times = Vector{FT}()
    data = Vector{FT}()
    for line in eachline(io)
        startswith(line, '\n') && continue
        line == "" && continue
        parsestart || begin
            parsestart |= startswith(line, 's')
            continue
        end
        t, abs = parseabsline(line, assert_)
        push!(times, t); push!(data, abs)
    end
    return AR{FT}(times, data)
end

function parseabsline(line, assert_=false)
    parsed = line |> split .|> q->parse(FT, q)
    if assert_
        @assert length(parsed) == 2 "only two elements on data lines"
    end

    t, abs = parsed
    return t, abs
end

Base.getindex(ar::AR, I::UnitRange{Int}) = ARView(ar, I.start, I.stop)

"""
    getARView(ar::AR, tstart, tstop)

Get a view of the `AbsorbanceRaw` object from a starting time to a stop time.
This can be useful for limiting analysis to a subset of the original data.
"""
function getARView(ar::AR, tstart::Real, tstop::Real)
    # get indices from times for the ar argument
    istart = time2index(ar, tstart)
    istop = time2index(ar, tstop)
    return ar[istart:istop]
end

const ARData = Union{AR, ARView}

