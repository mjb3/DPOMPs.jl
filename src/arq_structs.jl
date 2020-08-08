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
    sample_resolution::Int64
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
