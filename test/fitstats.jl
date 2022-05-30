# test/fitstats.jl --- series set and set group fitting and statistics

## just see if the functions don't error

@testset "Series set fitting" for setname in setnames
    ## just see it it works
    # construct a series set... (defined in utils.jl)
    # ... fit using internal functions
    results = fit(SeriesSetResults, getseriesset(setname))

    # trivial tests
    @test results isa SeriesSetResults
    @test length(results.concentrations) == length(results.initrates)
end

@testset "Set group fitting" begin
    
end
