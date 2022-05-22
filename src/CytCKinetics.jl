# By Bon Leif AMALLA (B4, 2022)

# Written for basic activity measurements. In using this initial version, there will
# be certain restrictions on the format of the input files and their filenames.
# More details can be found in the README file.

## PRELIM MODULE FOR TYPE DEFINITIONS AND READING ABSTRACTIONS

module CytCKinetics

include("absorbance_types.jl")
include("manip.jl")
include("order_fit.jl")
include("menten_kinetics.jl")

# series set calculations and manipulations
include("setmanips/basics.jl")
# include("setmanips/setmerge.jl")
include("setmanips/statistics.jl")

end # module

