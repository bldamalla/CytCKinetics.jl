# scripts for michaelis menten kinetics stuff

# the idea is to define a struct to which you store fit results and
# use that data for analysis.

# this script file defines an input dependent fitting function which can be
# modified to the format used

using CytCKinetics

"""
    quickfit(seriesdata; root)

Do a quick fit of the linear part of the series data `seriesdata`. This
argument is a `Dict` containing some specific keywords used in the configuration
file to be read.
"""
function quickfit(series; root)
    # open the file containing the data. Necssary to supply root folder
    data = open(joinpath(root, series["filename"])) do io
        read(io, AbsorbanceRaw)
    end

    viewstart, viewstop = series["viewrange"]
    tview = getARView(data, viewstart, viewstop)

    blanktime = series["blankat"] |> q->time2index(data, q)
    blanked = tview - data.values[blanktime]
    lned = log(blanked)

    # find appropriate N and corresponding r2 for fitting
    N, r2fitted = thresholdfit(lned; r2thresh=0.96, minthresh=20)

    # actual fitting
    k, _ = linfit(lned.times[1:N], lned.values[1:N])

    return FitResults{CytCKinetics.FT}(k, N, r2fitted)
end

"""
    mentenaxes(series; root)

Return the concentration and associated initial reaction rates given series
information.
"""
function mentenaxes(seriesdata; root)
    results = quickfit(seriesdata; root=root)

    conc = seriesdata["concentration"]
    # return the series concnetration and the fitted initial reaction rate
    return (conc, conc * -results.slope)
end

function ribbonfuncbase(x, p, covmatrix)
    ∇ = CytCKinetics.menten_jac(x, p)
    (∇ * covmatrix * transpose(∇))[] |> sqrt
end

@info "Loaded kinetics commons..."
