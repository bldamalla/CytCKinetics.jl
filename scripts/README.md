# Scripts README

The purposes of these scripts are to:
1. Show examples of module usage; and
2. For the main author to skip most REPL usage in analyses.

Lately, I am not motivated to write a documentation for a project that may
not even be shared. However, since transparency is essential I decided to
document how the module is used. Writing scripts also makes repeating
processes without bloating the main module with the author-specific and
spectrometer-specific formats.

The folder contains three files:
1. `setresults.jl` for calculating and fitting reaction rates
2. `setplotting.jl` for plotting reaction rates against concentrations
3. `sysimg.jl` (may be removed) for creating `CairoMakie` system images for faster
load time.

Files are run through VSCode during the development of the package. So far, there
are no plans to make the files accept command line arguments for the sets to be
analyzed since it's working anyway.
