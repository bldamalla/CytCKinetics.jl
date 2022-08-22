# CytCKinetics

A small Julia project containing utility functions for reading and calculation of kinetics related
constants. This package was specifically made for research use.

The package is a bare minimum of functional types that were thought to make analyses easier to
express. These will be expected to change in the near future once tests on further data are
possible. Other capabilities may also be added.

Existing capabilities include:
- Reading specific output format from a Hitachi UV-Vis spectrometer
- Utility types for handling absorbance data (including views and manipulated types)
- Simple data manipulations by extending methods from `Base`, to the whole data (like vectorized functions)
- Predefined model functions and Jacobians for Michaelis-Menten analyses using `LsqFit`
- Built-in methods for nonlinear fitting

Extending from these is not very difficult. A notebook/docs repository will be provided to show how
these extensions and analyses are done. They include:
- Checking reaction order
- Calculation of confidence intervals and hypothesis testing
- Extraction of some reaction steady-state properties

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

### Kinetic assay setup

Kinetic assays done in the lab usually use a fixed amount of enzyme while varying the initial
substrate concentration. The decay of absorbance at a certain wavelength is monitored and served as
a stand-in for the remaining substrate concentration. One assay uses a _set_ of initial substrate
concentrations; the absorbance time series of the measurements form a **series set**.

Michaelis--Menten kinetic parameters can be obtained from series sets. Series sets using the same
protein samples (done on the same day) can be combined to from **set groups** to get more degrees of
freedom (points). This _can_ give narrower confidence intervals for kinetic parameters.

Two modes of analyses coming from these points of view can be made and are handled interestingly
(yet still quite sloppily) by the package.

## Examples in notebooks and documentation

Using script files used for actual data analysis as examples for package has been discontinued. This
may increase the size of commits a lot. Instead, a **separate** repository containing both notebooks
and actual data used for analysis will be made not far into the future. As the repo will have data
I have obtained from experiments, it will use a different license from the one this repo uses.

Included with the notebooks is the documentation on how to use the package. Docstrings are already
being slowly added into the important constructs, so it will be easy to access them using the `help`
prompt in REPL. LSP implementations are also smart enough to cover these.

## Questions and contact

Send an email (see `Project.toml`) for questions regarding use or just about anything.

