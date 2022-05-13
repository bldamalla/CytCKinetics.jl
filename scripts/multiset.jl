# reading and calculating basic information for multiple sets
# these sets may come from the same strain or from different strains
# more documentation can be found in `scripts/README.md`

# there are no restrictions on the number of arguments
# but there are restrictions on how they are written

## parse the commandline arguments -- return an array of tuples
@assert length(ARGS) > 1 "use single series set analysis instead"
include("commons.jl")

# parse here
strainset_pairs = parse_setargument.(ARGS)

results = open(resultspath) do io
    TOML.parse(io)
end

for (strain, set) in strainset_pairs
    @assert haskey(results, strain) && haskey(results[strain], set) "Results for $(strain):$(set) not found in config file"
end

seriesdataset = map(strainset_pairs) do (strain, set)
    setresults = results[strain][set]
    fparams = setresults["fitparams"] isa Vector ? begin
        Vmax, kM = setresults["fitparams"]
        (Vmax, kM)
    end : nothing
    SeriesData(
        setresults["concentrations"],
        setresults["initialrates"],
        fparams,
        string(strain),
        string(set)
    )
end

@info "Set fitting results obtained... Setting up CairoMakie for plotting"
using CairoMakie, Colors

# set up empty plot and then add the information you need
fig = Figure()
ax = Axis(fig[1,1];
    xlabel="Initial [Cyt c (red)] (μM)",
    ylabel="Initial reaction rate (nM/s)",
    title="Reaction rates vs. concentration"
)

## TODO: make a recipe for this in commons I guess
## extend to more colors manually when needed
colors = range(Lab(60,-75,-75), stop=Lab(70,75,75), length=length(seriesdataset))
for (color, seriesdata) in zip(colors, seriesdataset)
    (; concentrations, initialrates) = seriesdata
    (; strain, setname) = seriesdata
    labelname = string(strain, ":", setname)
    scatter!(ax, concentrations, initialrates; label=labelname, color=color)

    mparams = seriesdata.mentenparams
    if !isnothing(mparams)
        xrange = 0:0.01:concentrations[end]
        Vmax, kM = mparams
        yrange = CytCKinetics.menten(xrange, [Vmax, kM])
        lines!(ax, xrange, yrange; color=color)
    end
end

Legend(fig[1,2], ax, framevisible=false)

savepath = joinpath(@__DIR__, "../plotdump", "multis", "tmp.pdf")
save(savepath, fig)
@info "Successfully produced simple graph in $(savepath)."
