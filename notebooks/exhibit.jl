### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# ╔═╡ 9da278f2-e8be-11ec-1ebb-c9e52e8b3a42
begin
	import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	Pkg.instantiate()
end;

# ╔═╡ 466b8b95-cf69-45e7-ba25-d9c5ac143b40
begin
	import TOML
	using CytCKinetics
end

# ╔═╡ 3285200d-46d5-4ccf-863d-bd5ee191377d
using CairoMakie

# ╔═╡ ac43f871-01f4-449f-92e8-10a14b81d8db
md"""
# Assay results and related notes

This notebook describes the steps taken using the module to
perform an analysis of an activity assay. This uses my data folder including
the configuration file format and the data format of the spectrometer.

Serialization and deserialization of results from the `results.toml` file
should not be included in the main module, but it has been included there
instead for convenience. For the purposes of this notebook, those functions
will not be used.
"""

# ╔═╡ ae178c97-d512-4340-b0ac-7f3292a53147
md"""
## Setting up

In the data folder, I have included a configuration file containing the
information the module uses for analyses. We first declare the "sets" we will
analyze. The format is `strain:set`.

This can be split into two parts using the below defined `split_arg` function.
As for the reason why the set format is such, it is a remnant of the previous
approach when I used command line arguments before a massive code refactor.
That approach may never be used anymore. The set format is still convenient, 
so that is still used.
"""

# ╔═╡ 7a01c188-6cc9-4d5d-95a7-47aef2e7d101
const dataloc = joinpath(@__DIR__, "..", "data")

# ╔═╡ bf7115ce-e06c-421f-bee3-c08942f16f00
setlist = ["WTC:Set1", "WTC:Set2"]

# ╔═╡ 08690860-9fe9-4dd9-a8e5-05c0ed5ec73d
split_arg(arg) = split(arg, ':') .|> string

# ╔═╡ fed398a4-be49-4705-b26f-babacc17c284
"""
	getsetmeta(strain, set)

A convenience function to extract the metadata for manipulating a set
of absorbance time series. Also returns the name of the folder where to
find absorbance data for the strain.
"""
function getsetmeta(strain, set)
	strainmeta = open(joinpath(@__DIR__, "../data", "$strain.toml")) do io
		TOML.parse(io)
	end
	root = strainmeta["rootfolder"]

	return strainmeta[set], root
end

# ╔═╡ e322b6e7-bf95-49a0-84ec-529deadff890
"""
	seriesfromdict(seriesmeta, root, sub)

A convenience function to create `AbsorbanceSeries` objects from the
extracted dictionary containing series information.
"""
function seriesfromdict(seriesmeta, root, sub)
	ardata = open(joinpath(dataloc, root, sub, seriesmeta["filename"])) do io
		read(io, AbsorbanceRaw)
	end

	framestart, framestop = seriesmeta["viewrange"]
	blanktime = seriesmeta["blankat"]

	return AbsorbanceSeries(ardata, framestart, framestop, blanktime)
end

# ╔═╡ b612d0ee-31ec-420e-9551-cc6c4184e93e
md"""
## Data extraction and fitting

We can finally extract data from a configuration file. Julia has a `TOML` module
in its standard library, so we can depend on that. We also need to parse input
from the spectrometer. The parser is included in the `CytCKinetics` module.

Absorbance time series are enclosed in `AbsorbanceSeries` objects containing
the raw absorbance values, the viewing frame from which to measure the reaction
rate and the perceived endpoint of the reaction. These objects are then collected
into `SeriesSet`s which are the basic units of analyses.

Regression can be performed on the `SeriesSet`s. Take note that fitting is an
opt-in operation (not enabled by default), so we must explicitly say that we
want the data set to be fit with the Michaelis-Menten model.
"""

# ╔═╡ b08d52c7-6268-4bf2-bd17-842f97f9074e
metadata = map(setlist) do strset
	strain, set = split_arg(strset)
	return getsetmeta(strain, set)
end

# ╔═╡ d857853e-014f-41f2-bf91-2f4f388d8c3f
seriessets = map(metadata) do setmeta
	meta, root = setmeta # root is rootfolder for strain data
	series_ = map(meta["Series"]) do seriesmeta
		seriesfromdict(seriesmeta, root, meta["subfolder"])
	end
	concentrations_ = map(meta["Series"]) do seriesmeta
		seriesmeta["concentration"]
	end
	menten_ = get(meta, "menten", false) # there has to be a menten=true setting

	return SeriesSet(series_, concentrations_, menten_, Dict{String,Any}())
end

# ╔═╡ a07280a8-76b3-44ab-9680-8a00346e8006
md"""
With the data converted to `SeriesSet`s, we can now proceed to fitting. The
`CytCKinetics` module exports the `SeriesSetResults` struct and the `fit`
function.

Fitting results are enclosed in `SeriesSetResults` objects. The `fit` function
was written in a way that we have to explicitly put in the type we want to
cast the results into. This is *not needed* but it helps with code clarity for
myself, hopefully also for other users.

The `fit` function accepts one optional position argument: `fitstart` which is
an estimate `[Vmax,Km]` for the kinetic parameters ``V_{max}`` and ``K_M``. It
also accepts two keyword arguments:
+ `r2thresh` = sets the threshold for the "first-order-ness" of the reaction rate for estimating the rate itself;
+ `minthresh` = sets the minimum number of time points needed for the reaction rate calculation.

For calculation details, consult the code (`thresholdfit` in `src/order_fit.jl`) or try to ask me directly.
"""

# ╔═╡ 28402550-1981-44cd-ae4d-8b80d782fbdf
fitresults = fit.(SeriesSetResults, seriessets; r2thresh=0.97, minthresh=20)

# ╔═╡ ecd16aa6-8729-4577-8383-c19df7da8669
md"""
## Plotting

I use Beamer to collect figures for printing, so the Makie backend installed is 
`CairoMakie`. Of course, any backend should work. Using Plots is also okay, but the 
plotting routine may be substantially different. There are plans to include plotting 
recipes in `Makie`, so until then there is no recommended plotting module.

In this notebook, information from all sets in `setlist` above are collected into one 
figure. This can be useful for comparison.
"""

# ╔═╡ 4b4a13f7-097b-4cfe-9432-2414a59cceec
begin
	f = Figure();
	ax = Axis(f[1,1],
		xlabel="Initial [Cyt c (red)] (μM)",
		ylabel="Initial reaction rate (nM/s)",
		title="Comparison of reaction rate curves across sets"
	);
end;

# ╔═╡ cadb3885-105f-452f-8d54-e3591ba9e712
for (i, (strset, fitresult)) in zip(setlist, fitresults) |> enumerate
	scatter!(ax, fitresult.concentrations, fitresult.initrates; label=strset, color=Cycled(i))

	if !isnothing(fitresult.fitparams)
		(; fitparams) = fitresult
		xvalues = 0.0:0.01:7.5
		yvalues = CytCKinetics.menten(xvalues, fitparams)
		lines!(ax, xvalues, yvalues; color=Cycled(i))
	end
end

# ╔═╡ a205fd8e-eebf-4538-a182-4caba6a222b7
leg = Legend(f[1,2], ax; framevisible=false);

# ╔═╡ 0d39fa38-1fad-41f5-945a-481f29e6347c
f

# ╔═╡ c14f6403-4d8e-4406-b64a-f93672d6c5a5
fparams = map(r->r.fitparams, fitresults) # Vmax and Km for the respective sets

# ╔═╡ e24580e5-ace6-4902-b0cb-3944e0ec88ad
md"""
The fitting parameters ``(V_{max}, K_M)`` for the curves are (from first to last):
+ ($(fparams[1][1]), $(fparams[1][2]))
+ ($(fparams[2][1]), $(fparams[2][2]))
"""

# ╔═╡ 32c2f422-cd1d-445a-ad67-0561bbb077ad
md"""
For context, the data used for generating the plots above were taken using **the
same** cytochrome *c* and oxidase samples. Ideally, the curves should be the same,
but since there may be a lot of errors between measurements, the curves are a bit
far off from each other.

In my opinion, they **should not be compared** with each other because they are
intended to be **replicates** of each other. Instead, the probably better approach
would be to combine the two sets into one "set group".

Inferential statistical methods are likely not applicable for sample sizes this small.
A good starting size may be 15 points.
"""

# ╔═╡ 2b4c3832-138a-4884-a8a2-09c36d05dc48
md"""
## Merging and comparison

Judging from the graphs, merging the two sets into a group makes the fitted curve
to be somewhere in the middle of the two curves. Merging can be done through the
internal `setmerge` function.
"""

# ╔═╡ 803cb185-ba9f-4f74-8929-7310e69d7162
mergedgroup = CytCKinetics.setmerge(SetGroup, seriessets...)

# ╔═╡ 701b1ba8-8341-4c1a-adeb-848ff7bc0383
md"""
Fitting is mostly identical. It should be kept in mind that `SetGroup`
objects are always subject to Michaelis-Menten fitting.
"""

# ╔═╡ 3126b51a-7468-483d-9bf7-94f7e1caebd1
groupfit = fit(SetGroupResults, mergedgroup)

# ╔═╡ 72e1370a-5e44-4477-b32b-8a30e0cf811e
md"""
We also create a new figure with all three curves in it for comparison.
"""

# ╔═╡ c8d46a43-5ea2-47b3-9146-8d83f5f73f77
begin
	g = Figure();
	ax2 = Axis(g[1,1],
		xlabel="Initial [Cyt c (red)] (μM)",
		ylabel="Initial reaction rate (nM/s)",
		title="Comparison of reaction rate curves across sets"
	);
end;

# ╔═╡ e6ae60bb-6a39-4140-bf34-bde50c7141f0
for (i, (strset, fitresult)) in zip(setlist, fitresults) |> enumerate
	scatter!(ax2, fitresult.concentrations, fitresult.initrates; label=strset, color=Cycled(i))

	if !isnothing(fitresult.fitparams)
		(; fitparams) = fitresult
		xvalues = 0.0:0.01:7.5
		yvalues = CytCKinetics.menten(xvalues, fitparams)
		lines!(ax2, xvalues, yvalues; color=Cycled(i))
	end
end

# ╔═╡ 22d47386-2bb9-45b0-b840-2398e1ca2983
let xvalues = 0.0:0.01:7.5
	yvalues=  CytCKinetics.menten(xvalues, groupfit.fitparams)
	lines!(ax2, xvalues, yvalues; color=Cycled(3), label="Merged group")
end;

# ╔═╡ 7086acbe-7e97-4efb-bf16-c8a46c45d3ec
leg2 = Legend(g[1,2], ax2; framevisible=false)

# ╔═╡ 0f33197c-721d-4b7e-bc51-e759c5175cce
g

# ╔═╡ 1c92b071-7a5b-454c-b3ae-c2aea37385d6
md"""
As previously inferred, the reaction rate curve for the merged group lies in
the *middle* of the respective curves. Moreover, it is higher than that of the first
set because there are contributions from the second set with slightly higher rates.
While this is the case, likely since most of the degrees of freedom come from the 
first set, the new curve *appears* closer to it than to the second set.
"""

# ╔═╡ 47f7354b-d16f-4420-a78d-8a94f2019ce5
md"""
## More statistics from `SetGroup`s

Since rate measurements contain errors, it may be essential to study these.
Two ways to do these are to:
1. Draw the confidence bands for the regression curve ("Delta method");
2. Draw the confidence ellipsoids for the fitted parameters (approximation)

Truly, concentrations also have measurement errors. With the current methods, it
may be hard to study those. The best that can be done is to perform the dilutions
as carefully as possible.

For the confidence intervals, we need the `Distributions` module.
"""

# ╔═╡ 43244bc1-5718-4691-9196-aa8f2112494c
import Distributions

# ╔═╡ 5759f262-2ea7-4ed6-9480-e1cf72749678
md"""
### Confidence bands for the regression curve

As it turns out it is possible to approximate the variance of the reaction rates
through the so-called [Delta method](https://en.wikipedia.org/wiki/Delta_method#Multivariate_delta_method).

The formula in the linked article is provided internally through the `fitstderror`
function. Usage will be explored in creating the confidence bands. We first create
a new figure containing only the points and the regression curve from the fitted
parameters.
"""

# ╔═╡ 9f9f4add-68be-4969-9399-35cd2d527dea
begin
	h = Figure()
	ax3 = Axis(h[1,1],
		xlabel="Initial [Cyt c (red)] (μM)",
		ylabel="Initial reaction rate (nM/s)",
		title="Merged group regression curve and its confidence band"
	);
end;

# ╔═╡ 2fa4c3ef-b6a2-4ee0-80f8-bc46c7a326d7
let xvalues = 0.0:0.01:7.5
	yvalues = CytCKinetics.menten(xvalues, groupfit.fitparams)
	scatter!(ax3, groupfit.concentrations, groupfit.initrates; color=Cycled(3))
	lines!(ax3, xvalues, yvalues; color=Cycled(3))
	α = 0.05

	dofs = length(groupfit.concentrations) - 2
	rel_coef = Distributions.quantile(Distributions.TDist(dofs), (1-α/2))

	σ_fse = map(xvalues) do x CytCKinetics.fitstderror(x, groupfit) * rel_coef end
	band!(ax3, xvalues, yvalues .- σ_fse, yvalues .+ σ_fse; color=(get(Makie.ColorSchemes.viridis, 0.75), 0.5))
end;

# ╔═╡ 536aa300-70fc-451b-8e5a-6c4292ba0f5f
h

# ╔═╡ e15285ae-3822-49bf-b5bb-ae55aa5e14c7
md"""
Above plotted is the regression curve for the merged group and its ``95\%``
confidence band. Its interpretation is that the *true* reaction rate curve 
lies in that region with ``95\%`` chance.
"""

# ╔═╡ 17803abe-9d42-4d04-aa52-75142d4bcf58
md"""
### Confidence regions for parameters

Calculation of the confidence regions is described within these lecture
[notes](https://stat.ethz.ch/~stahel/courses/cheming/nlreg10E.pdf). We begin creating the plots as follows:
"""

# ╔═╡ fdcb6753-d8af-46e8-9e95-8bef119de736
begin
	l = Figure()
	ax4 = Axis(l[1,1],
		title="Confidence region for kinetic parameters",
		xlabel="Vmax (nM/s)", ylabel="Km (mM)"
	)
end;

# ╔═╡ d22c2dbf-371a-4351-8a1a-2de82c72177f
let centers = groupfit.fitparams, widths = 3.5 .* groupfit.stderrors
	θ₁range = LinRange(centers[1]-widths[1], centers[1]+widths[1], 100)
	θ₂range = LinRange(centers[2]-widths[2], centers[2]+widths[2], 100)
	α = 0.05; α2 = 0.01
	
	dofp, dofn_p = (2, length(groupfit.concentrations)-2)
	q = Distributions.quantile(Distributions.FDist(dofp, dofn_p), 1-α)
	q2 = Distributions.quantile(Distributions.FDist(dofp, dofn_p), 1-α2)

	(; concentrations, initrates) = groupfit
	middle = CytCKinetics.fitssr(concentrations, initrates, centers)
	testfxn(θ₁, θ₂) = CytCKinetics.fitssr(concentrations, initrates, (θ₁, θ₂)) / middle
	θ₃vals = [testfxn(θ₁, θ₂) for θ₁ in θ₁range, θ₂ in θ₂range]

	scatter!(ax4, centers)
	contour!(ax4, θ₁range, θ₂range, θ₃vals; levels=[1 + dofp/dofn_p * q, 1 + dofp/dofn_p * q2])
end

# ╔═╡ 7c8505fd-64f1-4ff5-a1c3-11033fd7b94d
l

# ╔═╡ Cell order:
# ╟─9da278f2-e8be-11ec-1ebb-c9e52e8b3a42
# ╟─ac43f871-01f4-449f-92e8-10a14b81d8db
# ╟─ae178c97-d512-4340-b0ac-7f3292a53147
# ╠═7a01c188-6cc9-4d5d-95a7-47aef2e7d101
# ╠═bf7115ce-e06c-421f-bee3-c08942f16f00
# ╠═08690860-9fe9-4dd9-a8e5-05c0ed5ec73d
# ╟─fed398a4-be49-4705-b26f-babacc17c284
# ╟─e322b6e7-bf95-49a0-84ec-529deadff890
# ╟─b612d0ee-31ec-420e-9551-cc6c4184e93e
# ╠═466b8b95-cf69-45e7-ba25-d9c5ac143b40
# ╠═b08d52c7-6268-4bf2-bd17-842f97f9074e
# ╠═d857853e-014f-41f2-bf91-2f4f388d8c3f
# ╟─a07280a8-76b3-44ab-9680-8a00346e8006
# ╠═28402550-1981-44cd-ae4d-8b80d782fbdf
# ╟─ecd16aa6-8729-4577-8383-c19df7da8669
# ╠═3285200d-46d5-4ccf-863d-bd5ee191377d
# ╠═4b4a13f7-097b-4cfe-9432-2414a59cceec
# ╠═cadb3885-105f-452f-8d54-e3591ba9e712
# ╠═a205fd8e-eebf-4538-a182-4caba6a222b7
# ╟─0d39fa38-1fad-41f5-945a-481f29e6347c
# ╟─e24580e5-ace6-4902-b0cb-3944e0ec88ad
# ╟─c14f6403-4d8e-4406-b64a-f93672d6c5a5
# ╟─32c2f422-cd1d-445a-ad67-0561bbb077ad
# ╟─2b4c3832-138a-4884-a8a2-09c36d05dc48
# ╠═803cb185-ba9f-4f74-8929-7310e69d7162
# ╟─701b1ba8-8341-4c1a-adeb-848ff7bc0383
# ╠═3126b51a-7468-483d-9bf7-94f7e1caebd1
# ╟─72e1370a-5e44-4477-b32b-8a30e0cf811e
# ╠═c8d46a43-5ea2-47b3-9146-8d83f5f73f77
# ╠═e6ae60bb-6a39-4140-bf34-bde50c7141f0
# ╠═22d47386-2bb9-45b0-b840-2398e1ca2983
# ╠═7086acbe-7e97-4efb-bf16-c8a46c45d3ec
# ╟─0f33197c-721d-4b7e-bc51-e759c5175cce
# ╟─1c92b071-7a5b-454c-b3ae-c2aea37385d6
# ╟─47f7354b-d16f-4420-a78d-8a94f2019ce5
# ╠═43244bc1-5718-4691-9196-aa8f2112494c
# ╟─5759f262-2ea7-4ed6-9480-e1cf72749678
# ╠═9f9f4add-68be-4969-9399-35cd2d527dea
# ╠═2fa4c3ef-b6a2-4ee0-80f8-bc46c7a326d7
# ╟─536aa300-70fc-451b-8e5a-6c4292ba0f5f
# ╟─e15285ae-3822-49bf-b5bb-ae55aa5e14c7
# ╟─17803abe-9d42-4d04-aa52-75142d4bcf58
# ╠═fdcb6753-d8af-46e8-9e95-8bef119de736
# ╠═d22c2dbf-371a-4351-8a1a-2de82c72177f
# ╟─7c8505fd-64f1-4ff5-a1c3-11033fd7b94d
