## model comparison



"""
    run_model_comparison_analysis(model, obs_data; ... )

Run `n_runs` independent analyses for each `DPOMPModel` element in `models`, and compute [estimate] the **Bayesian model evidence** (BME.)

Returns an object of type `ModelComparisonResults`, which includes the mean and standard deviation of the estimates obtained.

**Parameters**
- `models`          -- An `Array` of `DPOMPModel`.
- `obs_data`        -- An array of type `Observation`.

**Optional**
- `n_runs`          -- number of algorithm runs used to estimate the bme (and variance.)
- `algorithm`       -- `String` representing the inference method used for the analysis, e.g. "SMC2" for *SMC^2* (default); "MBPI" for *MBP-IBIS*; or "ARQ".

**Inference algorithm parameters**
- `np`              -- number of [outer, i.e. theta] particles used in IBIS procedures (doesn't apply to ARQ-MCMC.)
- `ess_rs_crit`     -- Effective sampling size (resampling )

**Example**
```@repl
# NB. define some models to compare, e.g. as m1, m2, etc
models = [m1, m2, m3]
results = run_model_comparison_analysis(models, y; n_runs = 10)
tabulate_results(results)   # show the results (optional)
```

"""
function run_model_comparison_analysis(models::Array{DPOMPModel, 1}, y::Array{Observation, 1}; n_runs = 3, algorithm = C_ALG_NM_SMC2
    , np::Int64 = algorithm == C_ALG_NM_SMC2 ? C_DF_SMC2_P : C_DF_MBPI_P
    , ess_rs_crit::Float64 = algorithm == C_ALG_NM_SMC2 ? C_DF_ESS_CRIT : C_DF_MBP_ESS_CRIT
    , npf::Int64 = C_DF_PF_P, n_props::Int64 = C_DF_MBPI_MUT)

    println("Running: ", n_runs, "-run ", length(models), "-model Bayesian evidence analysis (algorithm := ", algorithm, ")\n - please note: this may take a while...")
    ## set up inference algorithm
    function alg_smc2(mdl::HiddenMarkovModel)
        theta_init = rand(mdl.prior, np)
        return run_pibis(mdl, theta_init, ess_rs_crit, true, C_ACCEPTANCE_ALPHA, npf).bme[1]
    end
    function alg_mibis(mdl::HiddenMarkovModel)
        theta_init = rand(mdl.prior, np)
        return run_mbp_ibis(mdl, theta_init, ess_rs_crit, n_props, false, C_ACCEPTANCE_ALPHA)
    end
    # function alg_arq(mdl::HiddenMarkovModel)    # WIP
    #     run_arq_mcmc_analysis(mdl, theta_range::Array{Float64,2}; sample_resolution::Int64 = 30, sample_limit::Int64 = 3, n_chains::Int64 = 5, steps::Int64 = C_DF_MCMC_STEPS, burnin::Int64 = df_adapt_period(steps), tgt_ar::Float64 = 0.33, np::Int64 = 200, ess_crit = 0.3)
    # end

    if algorithm == C_ALG_NM_SMC2
        inf_alg = alg_smc2
    elseif (SubString(algorithm, 1, 4) == C_ALG_NM_MBPI || algorithm == "MIBIS")
        inf_alg = alg_mibis
    else
        println(" WARNING - algorithm unknown: ", algorithm, "\n - defaulting to SMC2")
        inf_alg = alg_smc2
    end
    ## run analysis
    bme = zeros(n_runs, length(models))
    mnames = String[]
    start_time = time_ns()
    for m in eachindex(models)
        mdl = get_private_model(models[m], y)
        println(" processing model m", m, ": ", mdl.model_name)
        for n in 1:n_runs
            print("  analysis ", n, " ")
            bme[n, m] = inf_alg(mdl)
        end
        push!(mnames, mdl.model_name)
    end
    ## process results
    output = ModelComparisonResults(mnames, bme, vec(Statistics.mean(bme; dims = 1)), vec(Statistics.std(bme; dims = 1)), n_runs, time_ns() - start_time)
    println("Analysis complete (total runtime := ", Int64(round(output.run_time / C_RT_UNITS)), "s)")
    return output
end
