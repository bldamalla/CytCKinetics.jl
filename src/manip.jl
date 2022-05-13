# By Bon Leif AMALLA (B4, 2022)

# Written for basic activity measurements. In using this initial version, there will
# be certain restrictions on the format of the input files and their filenames.
# More details can be found in the README file.

## MANIPULATING ABSORBANCE DATA

export ARManip
export ARAny

### probably necessary math
mathsymbs_binary = [ :+, :-, :*, :/ ]

mathsymbs_unary = [
    :exp, :log,
    :sin, :cos, :tan, :asin, :acos, :atan, :cot
]

struct ARManip{T} <: AbstractVector{Tuple{T,T}}
    times::Vector{T}
    values::Vector{T}

    function ARManip(times, values)
        @assert length(times) == length(values)
        T = promote_type(eltype(times), eltype(values))
        return new{T}(times, values)
    end
end
Base.size(arm::ARManip) = (length(arm.times),)

function Base.getindex(arm::ARManip, i::Int)
    @boundscheck checkbounds(arm.times, i)
    return @inbounds (arm.times[i], arm.values[i])
end

const ARAny = Union{AbsorbanceRaw, ARView, ARManip}

for symb in mathsymbs_unary
    @eval begin
        function Base.$(symb)(ard::ARAny)
            ARManip(ard.times, Base.$(symb).(ard.values))
        end
    end
end

for symb in mathsymbs_binary
    @eval begin
        function Base.$(symb)(ard::ARAny, val)
            ARManip(ard.times, Base.$(symb).(ard.values, val))
        end
    end
end
            

