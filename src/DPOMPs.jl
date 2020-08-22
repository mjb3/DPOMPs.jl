## Bayesian inference for DSS Partially Observed Process Models in Julia
# by Martin Burke (contact: martin.burke@bioss.ac.uk)
# this software is available on the GPL License v3.0

# HMM: mcmc, mbp, smc, ibis, etc
# ARQ: production ARQ MCMC code
"""
DPOMPs is a package for:

- Bayesian parameter inference, and
- Simulation of,
- Discrete-state-space Partially Observed Markov Processes, in Julia.
- It also includes automated tools for convergence diagnosis and analysis.
- Developed for Julia `1.0`.
- Author: Martin Burke (martin.burke@bioss.ac.uk)
- Date: 2020-08-06

"""
module DPOMPs

### resources
import Statistics
import Distributions
import LinearAlgebra
import DataFrames
import CSV
import PrettyTables
import UnicodePlots     # https://github.com/Evizero/UnicodePlots.jl
import StatsBase
import Random           # for bitrand() function in arq

### global constants
const C_DEBUG = false
const C_RT_UNITS = 1000000000
const C_LBL_BME = "-ln p(y)"
const C_ALG_NM_SMC2 = "SMC2"
const C_ALG_NM_MBPI = "MBPI"
const C_INF_DELTA = 0.0000000000000001
const MAX_TRAJ = 196000
const C_PR_SIGDIG = 3

## MCMC
const C_DF_MCMC_STEPS = 50000
const C_DF_MCMC_ADAPT = 0.2
const C_MCMC_ADAPT_INTERVALS = 10
const C_ACCEPTANCE_ALPHA = 1.002

## IBIS
const C_DF_MBPI_P = 10000
const C_DF_SMC2_P = 4000
const C_DF_PF_P = 200
const C_DF_ESS_CRIT = 0.3
const C_DF_MBPI_ESS_CRIT = 0.5
const C_DF_MBPI_MUT = 3

## ARQ
const C_DF_ARQ_SL = 1
const C_DF_ARQ_SR = 30
const C_DF_ARQ_MC = 5
const C_DF_ARQ_AR = 0.33
const C_DF_ARQ_JT = 0.0

df_adapt_period(steps::Int64) = Int64(floor(steps * C_DF_MCMC_ADAPT))
# default_arq_sl(grng::Array{Float64,2}) = size(grng, 1) > 2 ? 1 : 7

### public stuffs ###
export DPOMPModel, Particle, Event, Observation
export SimResults, ImportanceSample, RejectionSample, MCMCSample, ARQMCMCSample, ModelComparisonResults
export generate_model, generate_custom_model, partial_gaussian_obs_model
export gillespie_sim, run_mcmc_analysis, run_ibis_analysis, run_arq_mcmc_analysis, run_model_comparison_analysis
export plot_trajectory, plot_parameter_trace, plot_parameter_marginal, plot_parameter_heatmap, plot_model_evidence
export get_observations, tabulate_results, print_results
export run_custom_mcmc_analysis, generate_custom_particle

#### DSS-POMPs ####

## types ###
include("hmm_structs.jl")
import Base: isless
isless(a::Event, b::Event) = isless(a.time, b.time)
isless(a::Observation, b::Observation) = isless(a.time, b.time)

## common ###
include("hmm_cmn.jl")
## Gillespie simulation ###
include("hmm_sim.jl")
## MCMC ###
include("hmm_mbp.jl")
include("hmm_std.jl")
include("hmm_mcmc.jl")
## generalised HMM pf ###
include("hmm_pf_resample.jl")
include("hmm_particle_filter.jl")
## IBIS sampling
include("hmm_ibis.jl")
## model comparison
include("hmm_mcomp.jl")
## predefined models ###
include("hmm_examples.jl")
## utils (e.g. printing to file) ###
include("hmm_utils.jl")
## visualisation (trajectories) ###
include("hmm_visuals_uc.jl")

#### public interface ## ####

"""
    gillespie_sim(model, parameters; tmax = 100.0, num_obs = 5)

Run a Doob-Gillespie (DGA) simulation based on `model`.

Returns a SimResults type containing the trajectory and observations data, or an array of the same if `n_sims` > 1.

**Parameters**
- `model`       -- `DPOMPModel` (see [DCTMPs.jl models]@ref).
- `parameters`  -- model parameters.
**Optional**
- `tmax`        -- maximum time (default: 100.)
- `n_obs`       -- the number of observations to draw (default: 5.)
- `n_sims`      -- number of simulations to draw (default: 1.)

**Example**
```@repl
using DPOMPs
m = generate_model("SIR", [50, 1, 0])
x = DPOMPs.gillespie_sim(model, [0.005, 0.12])
println(DPOMPs.plot_trajectory(x))
```

"""
function gillespie_sim(model::DPOMPModel, parameters::Array{Float64, 1}; tmax::Float64 = 100.0, num_obs::Int64 = 5, n_sims::Int64 = 1)
    y = generate_observations(tmax, num_obs, length(model.initial_condition))
    mdl = get_private_model(model, y)
    if n_sims == 1
        print("Running: ", model.model_name, " DGA for θ := ", parameters)
        output = gillespie_sim(mdl, parameters, true)
        println(" - finished.")
        return output
    else
        print("Running: ", model.model_name, " DGA for θ := ", parameters, " x ", n_sims)
        output = Array{SimResult,1}(undef, n_sims)
        for i in eachindex(output)
            y = generate_observations(tmax, num_obs, length(model.initial_condition))
            output[i] = gillespie_sim(mdl, parameters, true)
        end
        println(" - finished.")
        return output
    end
end

# Otherwise the results of a single-chain analysis are returned, which include the Geweke test statistics computed for that analysis.
"""
    run_mcmc_analysis(model, obs_data; ... )

Run an `n_chains`-MCMC analysis using the designated algorithm (*MBP-MCMC* by default.)

The `initial_parameters` are sampled from the prior distribution unless otherwise specified by the user. A Gelman-Rubin convergence diagnostic is automatically carried out (for n_chains > 1) and included in the [multi-chain] analysis results.

**Parameters**
- `model`               -- `DPOMPModel` (see [DCTMPs.jl models]@ref).
- `obs_data`            -- `Observations` data.

**Optional**
- `n_chains`            -- number of Markov chains (default: 3.)
- `initial_parameters`  -- 2d array of initial model parameters. Each column vector correspondes to a single model parameter.
- `steps`               -- number of iterations.
- `adapt_period`        -- number of discarded samples.
- `mbp`                 -- model based proposals (MBP). Set `mbp = false` for standard proposals.
- `ppp`                 -- the proportion of parameter (vs. trajectory) proposals in Gibbs sampler. Default: 30%. NB. not required for MBP.
- `fin_adapt`           -- finite adaptive algorithm. The default is `false`, i.e. [fully] adaptive.
- `mvp`                 -- increase for a higher proportion of 'move' proposals. NB. not applicable if `MBP = true` (default: 2.)

**Example**
```@repl
y = x.observations                          # some simulated data
model = generate_model("SIR", [50, 1, 0])   # a model
results = run_mcmc_analysis(model, y; fin_adapt = true) # finite-adaptive MCMC
tabulate_results(results)                   # optionally, show the results
```

"""
function run_mcmc_analysis(model::DPOMPModel, obs_data::Array{Observation,1}; n_chains::Int64 = 3, initial_parameters = rand(model.prior, n_chains), steps::Int64 = C_DF_MCMC_STEPS, adapt_period::Int64 = Int64(floor(steps * C_DF_MCMC_ADAPT)), fin_adapt::Bool = false, mbp::Bool = true, ppp::Float64 = 0.3, mvp::Int64 = 3)
    mdl = get_private_model(model, obs_data)
    if mbp
        return run_mbp_mcmc(mdl, initial_parameters, steps, adapt_period, fin_adapt)
    else
        return run_std_mcmc(mdl, initial_parameters, steps, adapt_period, fin_adapt, ppp, mvp)
    end
end

## MBP IBIS algorithm
"""
    run_mbp_ibis_analysis(model, obs_data; ... )

Run an *MBP-IBIS* analysis based on `model`, and `obs_data` of type `Observations`.

**Parameters**
- `model`               -- `DPOMPModel` (see [DCTMPs.jl models]@ref).
- `obs_data`            -- `Observations` data.

**Optional**
- `np`                  -- number of particles (default = 4000.)
- `ess_rs_crit`         -- resampling criteria (default = 0.5.)
- `n_props`             -- MBP mutations per step (default = 3.)
- `ind_prop`            -- true for independent theta proposals (default = false.)
- `alpha`               -- user-defined, increase for lower acceptance rate targeted (default = 1.002.)

**Example**
```@repl
# NB. using 'y' and 'model' as above
results = run_mbp_ibis_analysis(model, y)# importance sample
tabulate_results(results)                # show the results
```

"""
function run_mbp_ibis_analysis(model::DPOMPModel, obs_data::Array{Observation,1}; np = C_DF_MBPI_P, ess_rs_crit = C_DF_MBPI_ESS_CRIT, n_props = C_DF_MBPI_MUT, ind_prop = false, alpha = C_ACCEPTANCE_ALPHA)
    mdl = get_private_model(model, obs_data)
    theta_init = rand(mdl.prior, np)
    # (model::HiddenMarkovModel, theta_init::Array{Float64, 2}, ess_rs_crit = C_DF_ESS_CRIT; n_props = 3, ind_prop = false, alpha = 1.002)
    return run_mbp_ibis(mdl, theta_init, ess_rs_crit, n_props, ind_prop, alpha)
end

#### SMC ####
"""
    run_smc2_analysis(model, obs_data; ... )

Run an *SMC^2* (i.e. particle filter IBIS) analysis based on `model` and `obs_data` of type `Observations`.

**Parameters**
- `model`               -- `DPOMPModel` (see [DCTMPs.jl models]@ref).
- `obs_data`            -- `Observations` data.

**Optional**
- `np`                  -- number of [outer, i.e. theta] particles (default = 2000.)
- `npf`                 -- number of [inner] particles (default = 200.)
- `ess_rs_crit`         -- resampling criteria (default = 0.5.)
- `ind_prop`            -- true for independent theta proposals (default = false.)
- `alpha`               -- user-defined, increase for lower acceptance rate targeted (default = 1.002.)

**Example**
```@repl
# NB. using 'y' and 'model' as above
results = run_smc2_analysis(model, y)   # importance sample
tabulate_results(results)               # show the results
```
"""
function run_smc2_analysis(model::DPOMPModel, obs_data::Array{Observation,1}; np = C_DF_SMC2_P, npf = C_DF_PF_P, ess_rs_crit = C_DF_ESS_CRIT, ind_prop = true, alpha = C_ACCEPTANCE_ALPHA)
    mdl = get_private_model(model, obs_data)
    theta_init = rand(mdl.prior, np)
    # run_pibis(model::HiddenMarkovModel, theta::Array{Float64, 2}, ess_rs_crit::Float64, ind_prop::Bool, alpha::Float64, np::Int64
    println("Running: ", np, "-particle SMC^2 analysis (model: ", model.model_name, ")")
    return run_pibis(mdl, theta_init, ess_rs_crit, ind_prop, alpha, npf)
end

## NEW MAIN INTERFACE
"""
    run_ibis_analysis(model, obs_data; ... )

Run an iterative-batch-importance-sampling (IBIS) analysis based on `model` and `obs_data` of type `Observations`.

The default algorithm is *SMC^2* (i.e. particle filter IBIS), the other option is *model-based-proposal IBIS* (use `algorithm = "MBPI"`.) Note that the default value of the optional parameters below is sometimes affected by choice of algorithm. However these are overridden when specified by the user.

**Parameters**
- `model`               -- `DPOMPModel` (see [DCTMPs.jl models]@ref).
- `obs_data`            -- `Observations` data.

**Optional parameters**
- `algorithm`           -- `String`, per above.
- `np`                  -- number of [outer, i.e. theta] particles.

**Algorithm options**
- `ess_rs_crit`         -- resampling criteria.
- `ind_prop`            -- true for independent theta proposals.
- `alpha`               -- *increase* to *lower* the targeted acceptance rate (default = 1.002.)
- `npf`                 -- number of [inner] particles for *SMC^2* (default = 200.)
- `n_props`             -- MBP mutations per step (default: 3.)

**Example**
```@repl
# NB. using 'y' and 'model' as above
results = run_ibis_analysis(model, y)   # importance sample
tabulate_results(results)               # show the results
```
"""
function run_ibis_analysis(model::DPOMPModel, obs_data::Array{Observation,1}; algorithm::String = C_ALG_NM_SMC2
    , np = algorithm == C_ALG_NM_SMC2 ? C_DF_SMC2_P : C_DF_MBPI_P, ind_prop = algorithm == C_ALG_NM_SMC2
    , ess_rs_crit = algorithm == C_ALG_NM_SMC2 ? C_DF_ESS_CRIT : C_DF_MBPI_ESS_CRIT
    , alpha = C_ACCEPTANCE_ALPHA, npf = C_DF_PF_P, n_props = C_DF_MBPI_MUT)

    mdl = get_private_model(model, obs_data)
    theta_init = rand(mdl.prior, np)
    if algorithm == C_ALG_NM_SMC2
        println("Running: ", np, "-particle SMC^2 analysis (model: ", model.model_name, ")")
        return run_pibis(mdl, theta_init, ess_rs_crit, ind_prop, alpha, npf)
    else
        return run_mbp_ibis(mdl, theta_init, ess_rs_crit, n_props, ind_prop, alpha)
    end
end


#### ARQ-MCMC ####

## constants
const C_ALG_STD = "ARQ"
const C_ALG_DAUG  = "DAQ"
# const C_ALG_AD  = "ADARQ"

## structs
include("arq_structs.jl")
## visualisation tools
include("arq_visualisation_uc.jl")
### algorithms
## common functions, macro
include("arq_alg_cmn.jl")
## standard ARQ MCMC algorithm
include("arq_alg_std.jl")
## delayed acceptance ARQ MCMC algorithm
# include("arq_alg_da.jl")
## data augmented ARQ MCMC algorithm
include("arq_alg_daug.jl")
## common functions, printing, etc
include("arq_utils.jl")
## main algorithms
include("arq_main.jl")

end # module
