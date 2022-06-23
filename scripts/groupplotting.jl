## scripts/groupplotting.jl --- plotting with group information
## statistical happiness and hardships

using CairoMakie
import CytCKinetics, Distributions
import TOML

# get the groups you are interested in
arglist = ["WTC:GroupA", "WTC:GroupB"]
split_arg(arg) = split(arg, ':') .|> string

const resultsmetaloc = joinpath(@__DIR__, "..", "data", "results.toml")

resultsmeta = open(resultsmetaloc) do io
    TOML.parse(io)
end

# get the groupsetresults
results = map(arglist) do arg
    strain, group = split_arg(arg)
    return resultsmeta[strain][group] |> q->CytCKinetics.deserialize(CytCKinetics.SetGroupResults, q)
end
@info "Obtained results... Proceeding to plotting routine"

# plot all of the objects in the arglist i guess

f = Figure()
ax = Axis(f[1,1],
    title = "Reaction rate curves and their corresponding confidence bands",
    xlabel = "Initial [Cyt c (red)] (μM)",
    ylabel = "Initial reaction rate (nM/s)"
)
ax2 = Axis(f[2,1],
    title = "Confidence regions for kinetic parameters",
    ylabel = "Km (μM)", xlabel = "kcat (s⁻¹)"
)

paired = zip(arglist, results)
N = length(arglist)

for (i, (name, result)) in paired |> enumerate
    (; concentrations, initrates) = result
    gcolor = get(Makie.ColorSchemes.viridis, (i)/(N+1))
    scatter!(ax, concentrations, initrates;
        label=name, color=gcolor, markersize=5)

    # declaration of distribution parameters
    α = 0.05; dofN = length(result.concentrations) - 2; dofp = 2

    # the fitting curve and its confidence region
    let fparams = result.fitparams
        xvalues = 0:0.01:7.5
        yvalues = CytCKinetics.menten(xvalues, fparams)
        lines!(ax, xvalues, yvalues; color=gcolor)

        # the confidence band
        # declare a significance value of 0.05 for all experiments
        q = Distributions.quantile(Distributions.TDist(dofN), 1-α/2)
        
        σ_fse = map(xvalues) do x CytCKinetics.fitstderror(x, result) * q end
        band!(ax, xvalues, yvalues .- σ_fse, yvalues .+ σ_fse;
            color=(gcolor, 0.5))
    end

    # the confidence regions
    let centers = result.fitparams, widths = 4 .* result.stderrors
        θ1range = LinRange(centers[1]-widths[1], centers[1]+widths[1], 150)
        θ2range = LinRange(centers[2]-widths[2], centers[2]+widths[2], 150)
        
        q = Distributions.quantile(Distributions.FDist(dofp, dofN), 1-α)
        middlessr = CytCKinetics.fitssr(concentrations, initrates, centers)
        testfunc(iter) = CytCKinetics.fitssr(concentrations, initrates, iter) / middlessr
        θ3vals = [testfunc((θ1, θ2)) for θ1 in θ1range, θ2 in θ2range]

        scatter!(ax2, centers; color=gcolor, markersize=9, label=name)
        contour!(ax2, θ1range, θ2range, θ3vals; levels=[1+dofp/dofN * q], color = [gcolor])
    end

end

# legends
leg1 = Legend(f[1,2], ax; framevisible=false)
leg2 = Legend(f[2,2], ax2; framevisible=false)

const plotdump = joinpath(@__DIR__, "..", "plotdump/groups/thing.pdf")
save(plotdump, f)
@info "Saved image at $plotdump"
