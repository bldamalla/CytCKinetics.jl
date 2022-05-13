# command line args can include the strain to be measured and the specific set
# there should only be one argument

@assert length(ARGS) == 1 "single series set plotting requires 1 argument"
include("commons.jl")
include("menten_kinetics.jl")

## this assumes that all of the needed results have already been
## precalculated and stored in <datapath>/results.toml

strain, set = parse_setargument(ARGS[1])
setname = string(strain, set)

results = open(resultspath) do io
    TOML.parse(io)
end

@assert haskey(results, strain) && haskey(results[strain], set) "Results not found in config file"

setresults = results[strain][set]

# this script will only give out the plots for the seriesset
# all of the needed things for plotting should have already been
# calculated. It is recommended to use the `makie_img.so` sysimage
# when running julia

@info "Obtained set results... Setting up CairoMakie for plotting"

using CairoMakie, Colors

scatter_color = Lab(60,-75,-75)
line_color = Lab(70,75,75)

concentrations = setresults["concentrations"]
initialrates = setresults["initialrates"]

f = Figure()
ax = Axis(f[1,1];
    xlabel="Initial [Cyt c (red)] (μM)",
    ylabel="Initial reaction rate (nM/s)",
    title="Reaction rates vs. concentration"
)
scatter!(ax, concentrations, initialrates; label=setname, color=scatter_color)

if setresults["fitparams"] isa Vector
    fparams = setresults["fitparams"] |> SVector{2}
    xrange = 0:0.01:concentrations[end]
    yrange = CytCKinetics.menten(xrange, fparams)

    lines!(ax, xrange, yrange; color=line_color)
end

Legend(f[1,2], ax, framevisible=false)

savepath = joinpath(@__DIR__, "../plotdump/singles/$(strain)$(set).pdf")
save(savepath, f)
@info "Successfully saved figure at $savepath"
