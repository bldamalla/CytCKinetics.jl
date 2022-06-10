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

The package, as of writing, supports "series set analyses" wherein reaction rates are
calculated from a set of time series absorbance measurements and fit against corresponding
starting substrate concentrations. Plans are to include "set group analyses" for cleaner
statistical treatment not limited to parameter confidence ellipsoids and hypothesis tests.

## Scripts usage

So far, there are three Julia files in the `scripts/` folder. Their purposes are described
in the folder README. These files are edited depending on the data to be analyzed. Files
are run using VSCode, either in REPL mode or as a process

If there is enough motivation to do so, they may be edited in future version in ways that
allow them to accept command line arguments for flexibility.

## Questions and contact

Send an email (see `Project.toml`) for questions regarding use or just
about anything.
A documentation will be in progress once initial details are fleshed
out likely after the first semester (August 2022).
