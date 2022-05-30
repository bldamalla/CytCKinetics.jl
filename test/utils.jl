# test/utuls.jl --- utility functions based on some fixed format
# to read series set and set group information

function getseriesset(setname)
    # expect a string for setname

    metadata = open(configfile) do io
        TOML.parse(io)
    end

    setmeta = metadata[setname]
    seriesmeta = setmeta["Series"]
    seriescol = map(seriesmeta) do seriesdat
        serpath = joinpath(rawdatapath, setname, seriesdat["filename"])
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
    )

    concs = [serdat["concentration"] for serdat in seriesmeta]
    menten = setmeta["menten"]
    return SeriesSet(seriescol, concs, menten, setothermeta)
end
