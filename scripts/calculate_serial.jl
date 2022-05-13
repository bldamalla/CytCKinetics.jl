## calculation routine in serial (non-threaded)

using LsqFit

# this is just for the case that the user wants kinetic parameters
# recalculated for a series set... Defaults to false.
force_recalc = false

include("commons.jl")
include("menten_kinetics.jl")

## the script can be used for calculating and storing for multiple series sets
## though since code is written to be executed in serial, it may be better to
## use when there isn't a lot of new data (3 sets)

# expect the first command line argument to be about force recalc
if lowercase(ARGS[1]) == "true"
    sets = ARGS[2:end]
    if lowercase(ARGS[1]) == "true"
        force_recalc = true
    end
else
    sets = ARGS
end

# now look for the metadata
strainset_pairs = parse_setargument.(sets)
metadata_collection = [getmeta(strain) for (strain, _) in strainset_pairs]

metasetpairs = [(meta, set) for (meta, (_, set)) in zip(metadata_collection, strainset_pairs)]

rootfolders = [getsetroot(m, s) for (m, s) in metasetpairs]
seriessets = [getsetseries(m, s) for (m, s) in metasetpairs]

@info "Loaded series information metadata"

# now get the results... there is only one results file for everything
results = open(resultspath) do io
    TOML.parse(io)
end

# loop through the set names and see which has data or not
for (i, seriesset) in enumerate(seriessets)
    # check if the stuff can be found in results
    strain, set = strainset_pairs[i]
    straingrp = get(results, strain, nothing)
    if isnothing(straingrp) # then this is a completely new calculation for the strain
        results[strain] = Dict{String,Any}()
        straingrp = results[strain]
    end

    calculated = haskey(straingrp, set)

    # skip when already calculated in the past and there is no need to recalculate
    (calculated && !force_recalc) && continue
    @info "Setting up calculations for $(strain):$(set)"

    ## reaching here means you need to calculate kinetic parameters
    dataaxes = map(seriesset) do series
        mentenaxes(series; root=rootfolders[i])
    end
    concentrations = vcat([0.0], [concs for (concs, _) in dataaxes])
    initrates = vcat([0.0], [rates for (_, rates) in dataaxes]) .* 1000

    menten = let (meta, set) = metasetpairs[i]
        get(meta[set], "menten", false)
    end

    if menten
        @info "Michaelis-Menten fitting"

        fitmodel = curve_fit(CytCKinetics.menten_inplace,
                             CytCKinetics.menten_jac_inplace,
                             concentrations, initrates, fitstart;
                             inplace=true)
        fitparams = fitmodel.param
        stderrors = stderror(fitmodel)
        covvec = estimate_covar(fitmodel) |> vec
    else    # then there's no actual need for fitting; set to NaN since TOML doesnt
            # understand `nothing`
        fitparams = NaN
        stderrors = NaN
        covvec = NaN
    end

    ## now store as a dict and place the important information
    straingrp[set] = Dict(
        "concentrations" => concentrations,
        "initialrates" => initrates,
        "fitparams" => fitparams,
        "stderrors" => stderrors,
        "covvec" => covvec
    )
end

## create a fallback method jic there's only one
# string
cmp_(str) = identity(str)
cmp_(str::String, str2::String) = cmp(str, str2)

# rewrite the new results dict into the file
open(resultspath; write=true) do io
    TOML.print(io, results; sorted=true, by=cmp_)
end

@info "Printed results to $(resultspath)"
