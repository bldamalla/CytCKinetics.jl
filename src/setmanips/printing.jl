# By Bon Leif AMALLA (B4, 2022)

# Written for basic activity measurements. In using this initial version, there will
# be certain restrictions on the format of the input files and their filenames.
# More details can be found in the README file.

## CONVENIENCE CONVERSION OF DEFINED TYPES TO TOML COMPATIBLE ONES

function serialize(object::SeriesSetResults)
    retdict = Dict{String,Any}()

    retdict["concentrations"] = object.concentrations
    retdict["initrates"] = object.initrates

    # temporary cheating... if any of the last three are nothing
    # then fitting was not performed... so set toml value to NaN
    isnothing(object.stderrors) && begin
        retdict["fitparams"] = NaN
        retdict["stderrors"] = NaN
        retdict["covmatrix"] = NaN
        return retdict
    end

    # otherwise fitting was performed and then vectorize what needs
    # to be vectorized
    retdict["fitparams"] = vec(object.fitparams)
    retdict["stderrors"] = [object.stderrors...]
    retdict["covmatrix"] = vec(object.covmatrix)

    return retdict
end

function deserialize(::Type{SeriesSetResults}, dict)
    concentrations = dict["concentrations"]
    initrates = dict["initrates"]

    fitparams = dict["fitparams"] === NaN ? nothing : SVector{2}(dict["fitparams"])
    stderrors = dict["stderrors"] === NaN ? nothing : tuple(dict["stderrors"]...)
    covmatrix = dict["covmatrix"] === NaN ? nothing : SMatrix{2,2}(dict["covmatrix"])

    return SeriesSetResults(concentrations, initrates, fitparams, stderrors, covmatrix)
end
