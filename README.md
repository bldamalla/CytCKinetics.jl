# CytCKinetics

A small Julia project containing utility functions for reading and
calculation of kinetics related constants. This package was
specifically made for research use.

Data used for testing can be provided upon request.

The package is a bare minimum of functional types that were thought
to make analyses easier to express. These will be expected to change
in the near future once tests on further data are possible.
Other functionality may also be added.

Existing functionality include:
- Reading specific output format from a Hitachi UV-Vis spectrometer
- Utility types for handling absorbance data (including views and manipulated types)
- Simple data manipulations by extending methods from `Base`, to the whole data (like vectorized functions)
- Predefined model functions and Jacobians for Michaelis-Menten analyses using `LsqFit`

## Package usage

The module defines an `AbsorbanceRaw` type that contains time series
absorbance data obtained from a `.TXT` file input. The format of the
contents is likely specific to the software used by the spectrometer
and will not be shared. Though it is possible to recreate the
important parts of the data format based on the parser (quite simple!).

Raw data can be obtained as:
```julia
using CytCKinetics
ARdata = open("file.txt") do io
    read(io, AbsorbanceRaw)
end
```

Obtained data can be accessed elementwise similar to an `AbstractArray`.
Ranges can be accessed lazily giving `ARView` objects. "Broadcasted"
manipulations give rise to manipulated types `ARManip`. An example
of a manipulation is setting a blank from within the series:
```julia
startindex = time2index(ARdata, starttime)  # get index within data from time
stopindex = time2index(ARdata, stoptime)
blankindex = time2index(ARdata, blanktime)

viewdata = ARdata[startindex:stopindex]     # get data view
blankvalue = ARdata[blankindex][2]          # access blank absorbance
blanked = viewdata - blankvalue             # subtract the value from the entire view
                                            # creating a new ARManip
```

## Scripts usage

As of writing, scripts are tested to work against a computer running
on MacOS. The primary author also has a Windows machine for lab
related work and scripts are also being tested against that machine. So far, all
scripts work the same way (as expected, since Julia can be run cross-platform).
Useful scripts that utilize functions within the module are contained in the
`scripts/` folder.

So far, plans are to include two flavors of analyses:
1. Series set analyses; and
2. Set group analyses.
Set group analyses are the same as series set analyses, but with clearer statistical
treatment. The differences and actions of the scripts will be explained in the README within
the folder.

## Questions and contact

Send an email (see `Project.toml`) for questions regarding use or just
about anything.
A documentation will be in progress once initial details are fleshed
out likely after the first semester (August 2022).
