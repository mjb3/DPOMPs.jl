## for internal use (called by public functions)
function run_inner_mcmc_analysis(mdl::LikelihoodModel, da::Bool, steps::Int64, burnin::Int64, chains::Int64, tgt_ar::Float64)#, retain_samples::Bool
    ## for performance evaluation
    start_time = time_ns()
    ## designate inner MCMC function and results array
    mcmc_fn = da ? daarq_met_hastings! : arq_met_hastings!
    # is_mu_fn = da ? compute_da_is_mean : compute_is_mean
    # mcmc = Array{MCMCResults, 1}(undef, chains)
    n_theta = size(mdl.grid_range, 1)
    samples = zeros(n_theta, steps, chains)
    ## initialise importance sample
    grid = Dict() # {Array{Int64, 1}, GridPoint}
    # is_mu = zeros(n_theta)
    is_uc = 0.0
    ## run N chains
    for i in 1:chains
        ## retain_samples or clear IS
        # retain_samples || (grid = Dict())
        ## choose initial theta coords
        theta_init = rand(1:mdl.grid_resolution, n_theta)
        print(" initialising chain ", i, ": θ ~ ", round.(get_theta_val(mdl, theta_init, 0.0); sigdigits = C_PR_SIGDIG + 1))
        ## run inner MCMC using designated function
        # mcmc[i] = mcmc_fn(samples, i, grid, mdl, steps, burnin, theta_init, tgt_ar)
        mcmc = arq_met_hastings!(samples, i, grid, mdl, steps, burnin, theta_init, tgt_ar)
        println("- complete (calls to f(θ) := ", mcmc[1], "; AAR := ", round(mcmc[3] * 100, digits = 1), "%)")
    end
    ## compute scale reduction factor est.
    rejs = handle_rej_samples(samples, burnin)      # shared HMM functions
    sre = gelman_diagnostic(samples, burnin).sre
    ## get importance sample
    theta_w = collect_theta_weight(grid, n_theta)
    is_mu = zeros(n_theta)
    cv = zeros(length(is_mu),length(is_mu))
    # shared HMM fn:
    compute_is_mu_covar!(is_mu, cv, theta_w[1], theta_w[2]) #estimate_model_evidence(Statistics.mean(theta_w[2]))
    is_output = ImportanceSample(is_mu, cv, theta_w[1], theta_w[2], 0, [estimate_model_evidence(sum(theta_w[2]) / (mdl.grid_resolution ^ length(is_mu)))])
    ## return results
    output = ARQMCMCSample(is_output, rejs, burnin, mdl.grid_resolution, mdl.sample_limit, mdl.grid_range, sre, time_ns() - start_time)
    # output = ARQMCMCResults(da ? C_ALG_DA : C_ALG_STD, time_ns() - start_time, mdl.grid_resolution, mdl.sample_limit, mdl.jitter, steps, burnin, gmn[1], gmn[2], gmn[3], mdl.grid_range, grid, is_mu, estimate_model_evidence(is_uc), gmn[4], gmn[5], gmn[6], samples, mcmc)
    println("- finished in ", Int64(round(output.run_time / C_RT_UNITS)), "s. (Iμ = ", round.(is_output.mu; sigdigits = C_PR_SIGDIG), "; Rμ = ", round.(rejs.mu; sigdigits = C_PR_SIGDIG), "; BME = ", round.(output.imp_sample.bme[1]; sigdigits = C_PR_SIGDIG), ")")
    return output
end

## run standard ARQMCMC analysis
"""
    run_arq_mcmc_analysis(model::ARQModel; ... )

Run ARQMCMC analysis with `chains` Markov chains, where `n_chains > 1` the Gelman-Rubin convergence diagnostic is also run.

**Parameters**
- `model`               -- `ARQModel` (see docs.)
**Named parameters**
- `sample_resolution`   -- i.e. the length of each dimension in the importance sample.
- `sample_limit`        -- sample limit, should be increased when the variance of `model.pdf` is high (default: 1.)
- `n_chains`            -- number of Markov chains (default: 3.)
- `steps`               -- number of iterations.
- `burnin`              -- number of discarded samples.
- `tgt_ar`              -- acceptance rate (default: 0.33.)

"""
function run_arq_mcmc_analysis(model::ARQModel; sample_resolution::Int64 = 30, sample_limit::Int64 = 3, steps::Int64 = 50000, burnin::Int64 = Int64(floor(steps / 5)), n_chains::Int64 = 5, tgt_ar::Float64 = 0.33)#, retain_samples::Bool = true
    mdl = LikelihoodModel(model.pdf, model.parameter_range, sample_resolution, sample_limit, 0.0)
    println("Running: ARQMCMC analysis (", n_chains, " x " , steps, " steps):")
    return run_inner_mcmc_analysis(mdl, false, steps, burnin, n_chains, tgt_ar)
end

## - for direct access with internal model
function run_arq_mcmc_analysis(model::HiddenMarkovModel, theta_range::Array{Float64,2}; sample_resolution::Int64 = 30, sample_limit::Int64 = 3, n_chains::Int64 = 5, steps::Int64 = 50000, burnin::Int64 = Int64(floor(steps / 5)), tgt_ar::Float64 = 0.33, np::Int64 = 200, ess_crit = 0.3)
    mdl = ARQModel(get_log_pdf_fn(model, np; essc = ess_crit), theta_range)
    println(" ARQ model initialised: ", model.model_name)
    return run_arq_mcmc_analysis(mdl; sample_resolution = sample_resolution, sample_limit = sample_limit, n_chains = n_chains, steps = steps, burnin = burnin, tgt_ar = tgt_ar)
end

## - public interface for DPOMP.jl
"""
    run_arq_mcmc_analysis(model, ; ... )

Run ARQMCMC analysis with `chains` Markov chains, where `n_chains > 1` the Gelman-Rubin convergence diagnostic is also run.

**Parameters**
- `model`               -- `DPOMPModel` (see docs.)
**Named parameters**
- `sample_resolution`   -- i.e. the length of each dimension in the importance sample.
- `sample_limit`        -- sample limit, should be increased when the variance of `model.pdf` is high (default: 3.)
- `n_chains`            -- number of Markov chains (default: 3.)
- `steps`               -- number of iterations.
- `burnin`              -- number of discarded samples.
- `tgt_ar`              -- acceptance rate (default: 0.33.)
- `np`                  -- number of SMC particles in PF (default: 200.)
- `ess_crit`            -- acceptance rate (default: 0.33.)

"""
function run_arq_mcmc_analysis(model::DPOMPModel, obs_data::Array{Observation,1}, theta_range::Array{Float64,2}; theta_resolution::Int64 = 30, sample_limit::Int64 = 3, n_chains::Int64 = 5, steps::Int64 = 50000, burnin::Int64 = Int64(floor(steps / 5)), tgt_ar::Float64 = 0.33, np::Int64 = 200, ess_crit = 0.3)
    hmm = get_private_model(model, obs_data)
    return run_arq_mcmc_analysis(hmm, theta_range; sample_resolution = theta_resolution, sample_limit = sample_limit, n_chains = n_chains, steps = steps, burnin = burnin, tgt_ar = tgt_ar, np = np, ess_crit = ess_crit)
end

# ## run delayed acceptance ARQMCMC analysis
# """
#     run_daarq_mcmc_analysis(model, steps, adapt_period, chains::Int64 = 3)
#
# **Parameters**
# - `model`               -- `ARQModel` (see docs).
# - `sample_resolution`   -- i.e. the length of each dimension in the importance sample.
# - `steps`               -- number of iterations.
# - `burnin`              -- number of discarded samples.
# - `chains`              -- number of Markov chains (default: 3).
# **Named parameters**
# - `jitter`              -- add random noise to samples (0.0 to 0.5).
# - `da_limit`            -- delayed acceptance 'limit', i.e. threshold (default: 1).
# - `tgt_ar`              -- acceptance rate (default: 0.33).
# - `retain_samples`      -- for evaluation only, can be safely ignored (default: true).
#
# Run delayed acceptance ARQMCMC analysis with `chains` Markov chains. Where `chains > 1` the Gelman-Rubin convergence diagnostic is also run.
# """
# function run_daarq_mcmc_analysis(model::ARQModel, sample_resolution::Int64, steps::Int64, burnin::Int64, chains::Int64 = 3; jitter::Float64 = 0.25, da_limit::Int64 = 1, tgt_ar::Float64 = 0.33, retain_samples::Bool = true)
#     mdl = LikelihoodModel(model.pdf, model.parameter_range, sample_resolution, da_limit, jitter)
#     println("running DAARQ MCMC analysis (", chains, " x " , steps, " steps):")
#     return run_inner_mcmc_analysis(mdl, true, steps, burnin, chains, tgt_ar, retain_samples)
# end

# ## augmented data ARQ MCMC
# """
#     run_adarq_mcmc_analysis(model, steps, adapt_period, chains::Int64 = 3)
#
# **Parameters**
# - `model`               -- `ADARQModel` (see docs).
# - `sample_resolution`   -- i.e. the length of each dimension in the importance sample.
# - `steps`               -- number of iterations.
# - `burnin`              -- number of discarded samples.
# - `chains`              -- number of Markov chains (default: 3).
# **Named parameters**
# - `jitter`              -- add random noise to samples (0.0 to 0.5).
# - `da_limit`            -- delayed acceptance 'limit', i.e. threshold (default: steps).
# - `tgt_ar`              -- acceptance rate (default: 0.33).
# - `retain_samples`      -- for evaluation only, can be safely ignored (default: true).
#
# Run augmented data ARQ MCMC analysis with optional delayed acceptance (set `da_limit`). The Gelman-Rubin convergence diagnostic is run automatically.
# """
# function run_daq_mcmc_analysis(model::DAQModel, sample_resolution::Int64, steps::Int64, burnin::Int64, chains::Int64 = 3, jitter::Float64 = 0.25, da_limit::Int64 = steps, tgt_ar::Float64 = 0.33, retain_samples::Bool = true)
#     ## for performance evaluation
#     start_time = time_ns()
#     # MERGE SOME OF THIS STUFF? ***
#     mdl = LikelihoodModel(model.pdf, model.parameter_range, sample_resolution, da_limit, jitter)
#     println("running data augmented QMCMC analysis (", chains, " x " , steps, " steps):")
#     # return run_inner_mcmc_analysis(mdl, true, steps, burnin, chains, tgt_ar, retain_samples)
#     mcmc = Array{MCMCResults, 1}(undef, chains)
#     ## draw grid
#     # NEED TO ADD NEW STRUCTS
#     grid = Dict() # {Array{Int64, 1}, GridX}
#     is_mu = zeros(size(mdl.grid_range, 1))
#     ## run N chains
#     for i in eachindex(mcmc)
#         ## retain_samples && (grid = mcmc[i].grid)
#         retain_samples || (grid = Dict())
#         ## choose initial theta coords
#         theta_init = rand(1:mdl.grid_resolution, length(is_mu))    #length(DF_THETA)
#         theta_i = get_theta_val(mdl, theta_init)
#         print(" initialising chain ", i, ": θ = ", round.(theta_i; sigdigits = C_PR_SIGDIG + 1))
#         ## initialise x0
#         x0 = model.generate_x0(theta_i)
#         ## run (delayed acceptance) augmented data ARQ MCMC
#         mcmc[i] = adarq_met_hastings!(grid, mdl, x0, steps, burnin, theta_init, tgt_ar)   #, hsd, log_ir_y, prop_type
#         retain_samples || (is_mu .+= compute_da_is_mean(grid, length(is_mu)))
#         println(" - complete (AR: ", round(sum(mcmc[i].mc_accepted[burnin:steps]) / (steps - burnin) * 100, digits = 1), "%).")
#     end
#     ## compute scale reduction factor est.
#     gmn = gelman_diagnostic(mcmc, length(is_mu), steps - burnin)
#     ## compute importance sample mean
#     if retain_samples
#         is_mu .= compute_da_is_mean(grid, length(is_mu))[1]
#     else
#         is_mu ./= length(mcmc)
#     end
#     ## return results
#     output = ARQMCMCResults(C_ALG_DAUG, time_ns() - start_time, mdl.grid_resolution, mdl.sample_limit, mdl.jitter, steps, burnin, gmn[1], gmn[2], gmn[3], mdl.grid_range, grid, is_mu, 0.0, gmn[4], gmn[5], gmn[6], mcmc) #, prop_type
#     println("finished (sample μ = ", round.(output.mu; sigdigits = C_PR_SIGDIG), ").")
#     return output
# end
#
# function run_daq_mcmc_analysis(model::HiddenMarkovModel, theta_range::Array{Float64,2}, theta_resolution::Int64, steps::Int64, burnin::Int64, chains::Int64 = 3; jitter::Float64 = 0.25, da_limit::Int64 = steps, tgt_ar::Float64 = 0.33, retain_samples::Bool = true)
#     ## generate initial augmented data
#     function gen_x0(theta::Array{Float64,1})
#         x0::Particle = generate_x0(model, theta)
#         return AugDataResult(x0.log_like[2], x0)
#     end
#     ## evaluate pdf
#     function compute_density(x_i::AugDataResult, theta_f::Array{Float64,1})
#         ## make model based proposal
#         # current state
#         # pf = Discuit.ParameterProposal(theta_f, 1)  # C_DSCT_MODEL.prior_density(theta_f) - NOT USED CURRENTLY
#         mbp::Particle = model_based_proposal(model, theta_f, x_i.aug_data_var)
#         ## return as ARQMCMC type
#         return AugDataResult(mbp.log_like[2], mbp)
#     end
#     ## model
#     mdl = DAQModel(compute_density, gen_x0, theta_range)
#     return run_daq_mcmc_analysis(mdl, theta_resolution, steps, burnin, chains, jitter, da_limit, tgt_ar, retain_samples)
# end
