## model
#(see [Discuit.jl models](@ref) for further details)
"""
    ARQModel

**Fields**
- `pdf`         -- prior density function.
- `grid_range`  -- matrix representing the upper and lower limits of the parameter space.

A `struct` which includes the PDF (or estimate, or approximation of the target density - a function) and parameter space which specify the model.
"""
struct ARQModel{PFT<:Function}
    pdf::PFT
    parameter_range::Array{Float64, 2}
end

## augmented data model
"""
    DAQModel

**Fields**
- `pdf`         -- prior density function.
- `grid_range`  -- matrix representing the upper and lower limits of the parameter space.

Like `ARQModel` but for an augmented data model. The `pdf` has a signature like `pdf(xi::AugDataResult, theta::Array{Float64})` and must also return an `AugDataResult` (see the docs for further information).
"""
struct DAQModel{PFT<:Function, XFT<:Function}
    pdf::PFT
    generate_x0::XFT
    parameter_range::Array{Float64, 2}
end

## for internal use only
# NEED TO MERGE ARQModel?
struct LikelihoodModel{PFT<:Function}
    # model_name::String
    pdf::PFT
    grid_range::Array{Float64, 2}
    ## MOVE THESE TO algorithm functions?
    grid_resolution::Int64
    sample_limit::Int64
    jitter::Float64
end

## DA grid sample
struct GridSample
    sample::Array{Float64, 1}   # i.e. theta
    log_likelihood::Float64
end

## DA grid 'value'
struct GridSet
    anchor::Array{Float64, 1}   # i.e. theta REQ'D? ***
    samples::Array{GridSample, 1}
    # log_likelihood::Float64     # i.e weighted density
end

## DA grid request
struct DAGridRequest
    set::GridSet
    result::GridSample
    delayed::Bool
end

## grid point (internal)
# BENCHMARK with in place (NO UPDATING GRID ***)
struct GridPoint
    sample::Array{Float64, 1}   # i.e. theta
    log_likelihood::Float64
    visited::Int64  # i.e. # times PDF computed (max: N)
    sampled::Int64  # number of adapted samples
end

## results for a density estimate (internal)
struct GridRequest
    result::GridPoint
    process_run::Bool
end

## mcmc results data
# """
#     MCMCResults
#
# **Fields**
# - `process_count`   -- number of distinct process runs, i.e. f(x).
# - `run_time`        -- algorithm run time.
# - `grid`            -- sparse QMC structure.
# - `sample_idx`      -- two dimensional array of sample indices.
# - `samples`         -- two dimensional array of samples.
# - `mc_accepted`     -- proposal accepted array.
# - `mean`            -- sample mean.
# - `covar`           -- parameter covariance matrix.
# - `proposal_alg`    -- proposal algorithm.
# - `num_obs`         -- number of observations.
#
# The results of an individual ARQMCMC run. Includes samples; mean; covariance matrix; adaptation period; and results of the Geweke test of stationarity.
# """
# struct MCMCResults
#     process_count::Int64
#     sample_idx::Array{Int64, 2}
#     # samples::Array{Float64, 2}
#     mc_accepted::BitArray{1}
#     mc_fx::Array{Int16, 1}
#     mean::Array{Float64, 1}
#     covar::Array{Float64, 2}
#     # TO BE REMOVED (DEBUG)
#     mcf::Array{Float64, 2}  # TO BE REMOVED
#     mc_log_like::Array{Float64, 1}  # TO BE REMOVED
#     mc_time::Array{UInt64, 1}  # TO BE REMOVED
# end

## ARQMCMCResults, incl. Gelman Rubin test results
# """
#     ARQMCMCResults
#
# **Fields**
# - `sample_limit`    -- sample limit (i.e. max pdf calls per coordinate set).
# - `grid_resolution` -- grid resolution.
# - `adapt_period`    -- adaptation (i.e. 'burn in') period.
# - `mu`              -- between chain sample mean.
# - `sre`             -- scale reduction factor estimate (ADD XREF: see Gelman diagnostic).
# - `sre_ll`          -- scale reduction factor lower confidence interval.
# - `sre_ul`          -- scale reduction factor upper confidence interval.
# - `mcmc`            -- array of `MCMCResults`
#
# Results of a limited sample size quasi MCMC (ARQMCMC) analysis, including N `MCMCResults` variables; `mu`; and the scale reduction factor estimates (`sre`), i.e. Gelman Rubin convergence diagnostic.
# """
# struct ARQMCMCResults
#     algorithm::String
#     run_time::UInt64
#     grid_resolution::Int64
#     sample_limit::Int64
#     jitter::Float64
#     iterations::Int64
#     adapt_period::Int64
#     # parameters
#     # mu::Array{Float64, 1}
#     sdb::Array{Float64, 1}
#     sdw::Array{Float64, 1}
#     grid_range::Array{Float64, 2}
#     grid::Dict{Any, Any}
#     is_mu::Array{Float64, 1}
#     is_bme::Float64
#     # gelman diagnostic
#     sre::Array{Float64, 1}
#     sre_ll::Array{Float64, 1}
#     sre_ul::Array{Float64, 1}
#     # mcmc
#     samples::RejectionSample
#     mcmc::Array{MCMCResults,1}
# end

## autocorrelation results
"""
    AutocorrelationResults

**Fields**
- `lag`             -- autocorrelation lag.
- `autocorrelation` -- autocorrelation statistics.

Results of a call to `compute_autocorrelation`.
"""
struct AutocorrelationResults
    lag::Array{Int64,1}
    autocorrelation::Array{Float64,2}
end
