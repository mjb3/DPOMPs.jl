## event
"""
    Event

Requires no explanation.

**Fields**
- `time`        -- the time of the event.
- `event_type`  -- indexes the rate function and transition matrix.

"""
struct Event
    time::Float64
    event_type::Int64
end

## observation tuple
"""
    Observation

A single observation. Note that by default `val` has the same size as the model state space. However that is not necessary - it need only be compatible with the observation model.

**Fields**
- `time`        -- similar to `Event.time`, the time of the observation.
- `obs_id`      -- <1 if not a resampling step.
- `prop`        -- optional information for the observation model.
- `val`         -- the observation value.

"""
struct Observation
    time::Float64
    obs_id::Int64   # <1 if not a resampling step
    prop::Float64 #df: 1.0
    val::Array{Int64,1}
end

## a single realisation of the model
"""
    Particle

E.g. the main results of a simulation including the initial and final conditions, but not the full state trajectory.

**Fields**
- `theta`               -- e.g. simulation parameters.
- `initial_condition`   -- initial system state.
- `final_condition`     -- final system state.
- `trajectory`          -- the event history.
- `log_like`            -- trajectory log likelihood, mainly for internal use.

"""
struct Particle
    theta::Array{Float64,1}
    initial_condition::Array{Int64}
    final_condition::Array{Int64}
    trajectory::Array{Event,1}
    log_like::Array{Float64,1}  # log prior; full log like g(x); [latest] marginal g(x) / proposal likelihood (SMC / MCMC)
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
"""
    SimResults

The results of a simulation, including the full state trajectory.

**Fields**
- `model_name`      -- string, e,g, `"SIR"`.
- `particle`        -- the 'trajectory' variable, of type `Particle`.
- `population`      -- records the final system state.
- `observations`    -- simulated observations data (an `Array` of `Observation` types.)

"""
struct SimResults
    model_name::String
    particle::Particle
    population::Array{Array{Int64},1}
    observations::Array{Observation,1}
end

## public model
"""
    DPOMPModel

A `mutable struct` which represents a DSSCT model (see [Models](@ref) for further details).

**Fields**
- `model_name`          -- string, e,g, `"SIR"`.
- `rate_function`       -- event rate function.
- `initial_condition`   -- initial condition.
- `m_transition`        -- transition matrix.
- `obs_function         -- observation function, use this to add 'noise' to simulated observations.
- `obs_model`           -- observation model likelihood function.
- `prior`               -- prior [multivariate] Distributions.Distribution.
- `t0_index`            -- index of the parameter that represents the initial time. `0` if fixed at `0.0`.

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
"""
    RejectionSample

Essentially, the main results of an MCMC analysis, consisting of samples, mean, and covariance matrix.

**Fields**
- `samples`         -- three dimensional array of samples, e.g. parameter; iteration; Markov chain.
- `mu`              -- sample mean.
- `cv`              -- sample covariance matrix.

"""
struct RejectionSample
    theta::Array{Float64,3}         # dims: theta index; chain; sample
    mu::Array{Float64,1}
    cv::Array{Float64,2}
end

## IBIS sample
"""
    ImportanceSample

The results of an importance sampling analysis, such as iterative batch importance sampling algorithms.

**Fields**
- `mu`              -- weighted sample mean.
- `cv`              -- weighted covariance.
- `theta`           -- two dimensional array of samples, e.g. parameter; iteration.
- `weight`          -- sample weights.
- `run_time`        -- application run time.
- `bme`             -- Estimate (or approximation) of the Bayesian model evidence.

"""
struct ImportanceSample
    mu::Array{Float64,1}
    cv::Array{Float64,2}
    theta::Array{Float64,2}
    weight::Array{Float64,1}
    run_time::UInt64
    bme::Array{Float64,1}
end

## MBP MCMC, PMCMC
"""
    MCMCSample

The results of an MCMC analysis, mainly consisting of a `RejectionSample`.

**Fields**
- `samples`         -- samples of type `RejectionSample`.
- `adapt_period`    -- adaptation (i.e. 'burn in') period.
- `sre`             -- scale reduction factor estimate, i.e. Gelman diagnostic.
- `run_time`        -- application run time.

"""
struct MCMCSample
    samples::RejectionSample
    adapt_period::Int64
    sre::Array{Float64,2}
    run_time::UInt64
end

## ARQMCMC
"""
    ARQMCMCSample

The results of an ARQ MCMC analysis including the ImportanceSample and resampled RejectionSample.

The `sre` scale factor reduction estimates relate the rejection (re)samples to the underlying importance sample.

**Fields**
- `imp_sample`      -- main results, i.e. ImportanceSample.
- `samples`         -- resamples, of type RejectionSample.
- `adapt_period`    -- adaptation (i.e. 'burn in') period.
- `sample_resolution` -- number of distinct [possible] sample values along each dimension in the unit cube.
- `sample_limit`    -- maximum number of samples per theta tupple.
- `grid_range`      -- bounds of the parameter space.
- `sre`             -- scale reduction factor estimate, i.e. Gelman diagnostic. NB. *only valid for resamples*.
- `run_time`        -- application run time.

"""
struct ARQMCMCSample
    imp_sample::ImportanceSample
    samples::RejectionSample
    adapt_period::Int64
    sample_resolution::Int64
    sample_limit::Int64
    #     jitter::Float64
    grid_range::Array{Float64,2}
    sre::Array{Float64,2}
    run_time::UInt64 # redundant
end
