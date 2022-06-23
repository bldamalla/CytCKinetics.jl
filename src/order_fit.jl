# By Bon Leif AMALLA (B4, 2022)

# Written for basic activity measurements. In using this initial version, there will
# be certain restrictions on the format of the input files and their filenames.
# More details can be found in the README file.

## FITTING CURVE WITH PSEUDO FIRST ORDER KINETICS

using Statistics

export thresholdfit
export linfit

function thresholdfit(ar::ARAny, n=length(ar); r2thresh=0.975, minthresh=15)
    f = firstindex(ar)
    ts, abs = @inbounds ar.times[f:f-1+n], ar.values[f:f-1+n]

    # calculate the first r2 and then regress n while less than threshold
    curr_r2 = r2(ts, abs)
    n <= minthresh && return minthresh, begin
        ts, abs = @inbounds ar.times[f:f-1+minthresh], ar.values[f:f-1+minthresh]
        lastr2 = r2(ts, abs)
        @info "$minthresh"
        lastr2 < r2thresh && @warn "R² threshold ($r2thresh) not reached at $lastr2"
        lastr2
    end
    curr_r2 > r2thresh && return n, curr_r2

    # recalculate bounds for line calculation
    thresholdfit(ar, div(n, 8)*7; r2thresh=r2thresh, minthresh=minthresh)
end

function linfit(xs, ys)
    @assert length(xs) == length(ys) "data to be fit must have the same length"
    l = length(xs)

    matxs = hcat(xs, ones(l))
    m, b = matxs \ ys
    return m, b
end

rcoeff(xs, ys) = cov(xs, ys) / std(xs) / std(ys)
r2(xs, ys) = (rcoeff(xs, ys))^2

