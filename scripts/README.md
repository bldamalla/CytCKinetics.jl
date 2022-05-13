# Scripts README

As explained in the root README, this will contain explanations about what the
scripts do, especially statistical theory for series group data analyses.

There are two flavors for analysis of kinetics time series data from absorbance:
1) plain determination of reaction rates and fitting of Michaelis--Menten parameters from
sets of time series absorbance data (series set analyses); and
2) analyses of series sets belonging to a group and the respective statistics of fitting
parameters (set group analyses).

The second flavor is essentially a statistical extension of the first one. This allows
rigorous comparison between treatments, but one has to make sure (or assume) that sufficient
data has been obtained to properly conclude.

As of writing, the treatment for set group analyses has not been fleshed out yet to a certain
degree, and it is recommended to use series set analyses for the meantime. Qualitative and
non-rigorous quantitative comparisons can already be made between groups by looking at the 
fitting curves and comparing the confidence intervals for fitted parameters.

In general, scripts are general enough such that they can be used across operating systems that
can run Julia. However, priority support is for MacOS and Windows since the primary author only
has such machines. Also as to decrease the time for plotting (mainly Makie precompilation)
an additional script `sysimg.jl` has been included to create system images that use the
Makie routines deemed useful.

## Series set analyses

### Running scripts

Analyses are composed of the following steps:
1. Calculation of reaction rates from time-series absorbance data;
2. Calculation of Michaelis--Menten parameters;
3. Visualization routines (plotting with Makie)

The first (and, optionally, second steps) are done by running the following line in shell
(assuming that the working directory is the project root)
```sh
<juliaruntime> --project scripts/calculate_serial.jl <strain>:<set>...
```
The part `<juliaruntime>` means the path to the Julia executable. Normally, this can just be a
plain call to `julia`, but there are cases where the executable is not included in `PATH`.
The script then looks at the configuration files and looking for the `strain` and `set` to be
analyzed.

### About `results.toml` files

These files serve as intermediaries between calculation of reaction rates and visualization
routines, and are automatically generated. Rule of thumb is to not manually edit these. This
is fraud. Recalculations for series sets can also be forced.

## Set group analyses

<!-- TODO: write more about this once stat concepts are cleared out -->
