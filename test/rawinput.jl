# tests/rawinput.jl --- test if reading the input looks fine

import TOML

@testset "Reading series set input" begin
    strainmeta = open(joinpath(configfile)) do io
        TOML.parse(io)
    end
    stnm = setnames[1]

    @testset "Reading single series input" begin
        seriesmeta = strainmeta[stnm]["Series"][1]
        seriespath = joinpath(rawdatapath, stnm, seriesmeta["filename"])
        rawdata = open(seriespath) do io
            read(io, AbsorbanceRaw)
        end

        @test rawdata isa AbsorbanceRaw

        @test starttime(rawdata) == 0.0
        @test stoptime(rawdata) == 122.0
        @test timestep(rawdata) == 0.5
    
        ## define an absorbance series
        serstart, serstop = seriesmeta["viewrange"]
        blanktime = seriesmeta["blankat"]
        series = AbsorbanceSeries(rawdata, serstart, serstop, blanktime)

        @test series isa AbsorbanceSeries

        @test series.framestart == 39
        @test series.framestop == 57
        @test series.blanktime == 90.5
    end

    setmeta = strainmeta[stnm]
    seriesmetadata = setmeta["Series"]
    seriescol = map(seriesmetadata) do seriesdat
        serpath = joinpath(rawdatapath, stnm, seriesdat["filename"])
        rawdata = open(serpath) do io
            read(io, AbsorbanceRaw)
        end

        serstart, serstop = seriesdat["viewrange"]
        blanktime = seriesdat["blankat"]
        return AbsorbanceSeries(rawdata, serstart, serstop, blanktime)
    end

    excluded_keys = ["menten", "Series", "subfolder"]
    setothermeta = Dict(
        k=>setmeta[k] for k in keys(setmeta) if k ∉ excluded_keys
    ) # is there a more efficient way of doing this?
    concs = [serdat["concentration"] for serdat in seriesmetadata]
    menten = setmeta["menten"]
    serset = SeriesSet(seriescol, concs, menten, setothermeta)

    @test serset isa SeriesSet

    @test serset.menten === true
    @test length(serset.concentrations) == length(serset.series)
    @test serset.concentrations[begin] == 0.3
    @test serset.concentrations[end] == 7.5
    # probably assert this on creation instead of relying on user ??
end

@testset "Reading group metadata" begin
    
end
