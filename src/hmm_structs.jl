## event
struct Event
    time::Float64
    event_type::Int64
end

## observation tuple
struct Observation
    time::Float64
    obs_id::Int64   # <1 if not a resampling step
    prop::Float64 #df: 1.0
    val::Array{Int64,1}
end

## a single realisation of the model
struct Particle
    theta::Array{Float64,1}
    initial_condition::Array{Int64}
    final_condition::Array{Int64}
    trajectory::Array{Event,1}
    log_like::Array{Float64,1}  # log prior; full log like g(x); [latest] marginal g(x)
end
# - for dependent f/g
struct DFGParticle
    theta::Array{Float64,1}
    initial_condition::Array{Int64}
    final_condition::Array{Int64}
    trajectory::Array{Event,1}
    log_like::Array{Float64,1}  # prior, g(x)
    g_trans::Array{Int64,2}
end

## results of gillespie sim
struct SimResults
    model_name::String
    particle::Particle
    population::Array{Array{Int64},1}
    observations::Array{Observation,1}
end

## public model
"""
    DPOMPModel

**Fields**
- `model_name`          -- string, e,g, `"SIR"`.
- `initial_condition`   -- initial condition
- `rate_function`       -- event rate function.
- `m_transition`        -- transition matrix.
- `observation_function -- observation function, use this to add 'noise' to simulated observations.
- `prior_density`       -- prior density function.
- `observation_model`   -- observation model likelihood function.
- `t0_index`            -- index of the parameter that represents the initial time. `0` if fixed at `0.0`.

A `mutable struct` which represents a DSSCT model (see [Discuit.jl models](@ref) for further details).
"""
mutable struct DPOMPModel
    model_name::String                  # model name
    rate_function::Function             # computes event rates (in place)
    initial_condition::Array{Int64,1}   # sets (e.g. draw a sample from some known density) initial condition
    m_transition::Array{Int64,2}        # i.e adjusts the population according to event type
    obs_function::Function              # observation function (sim only) - TO BE REMOVED?
    obs_model::Function                 # observation model (log likelihood)
    prior::Distributions.Distribution   # prior distribution
    t0_index::Int64                     # == 0 if initial time known
end

## DAC private model
struct HiddenMarkovModel{RFT<:Function, ICT<:Function, TFT<:Function, OFT<:Function, OMT<:Function, PFT<:Distributions.Distribution}
    model_name::String                  # model name
    n_events::Int64                     # number of event types
    rate_function::RFT                  # computes event rates (in place)
    fn_initial_condition::ICT           # sets (e.g. draw a sample from some known density) initial condition
    fn_transition::TFT                  # i.e adjusts the population according to event type
    obs_function::OFT                   # observation function (sim only) - TO BE REMOVED
    obs_model::OMT                      # observation model (log likelihood)
    obs_data::Array{Observation,1}      # obs data
    # fn_log_prior::PFT               # prior density function (log likelihood)
    prior::PFT   # prior distribution
    t0_index::Int64                     # == 0 if initial time known
end

## generic rejection sample
struct RejectionSample
    theta::Array{Float64,3}         # dims: theta index; chain; sample
    mu::Array{Float64,1}
    cv::Array{Float64,2}
end

## IBIS sample
struct ImportanceSample
    mu::Array{Float64,1}
    cv::Array{Float64,2}
    theta::Array{Float64,2}
    weight::Array{Float64,1}
    run_time::UInt64
    bme::Array{Float64,1}
end

## MBP MCMC, PMCMC
struct MCMCSample
    samples::RejectionSample
    adapt_period::Int64
    sre::Array{Float64,2}
    run_time::UInt64
end

## ARQMCMC
# ADD process count? **
struct ARQMCMCSample
    imp_sample::ImportanceSample
    samples::RejectionSample
    adapt_period::Int64
    grid_resolution::Int64
    sample_limit::Int64
    #     jitter::Float64
    grid_range::Array{Float64,2}
    sre::Array{Float64,2}
    run_time::UInt64 # redundant
end
