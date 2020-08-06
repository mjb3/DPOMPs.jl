var documenterSearchIndex = {"docs":
[{"location":"manual/#DPOMPs.jl-manual","page":"DPOMPs.jl manual","title":"DPOMPs.jl manual","text":"","category":"section"},{"location":"manual/","page":"DPOMPs.jl manual","title":"DPOMPs.jl manual","text":"Pages = [\"manual.md\"]\nDepth = 3","category":"page"},{"location":"manual/#Types","page":"DPOMPs.jl manual","title":"Types","text":"","category":"section"},{"location":"manual/","page":"DPOMPs.jl manual","title":"DPOMPs.jl manual","text":"DPOMPModel\nSimResults","category":"page"},{"location":"manual/#DPOMPs.DPOMPModel","page":"DPOMPs.jl manual","title":"DPOMPs.DPOMPModel","text":"DPOMPModel\n\nFields\n\nmodel_name          – string, e,g, \"SIR\".\ninitial_condition   – initial condition\nrate_function       – event rate function.\nm_transition        – transition matrix.\n`observation_function – observation function, use this to add 'noise' to simulated observations.\nprior_density       – prior density function.\nobservation_model   – observation model likelihood function.\nt0_index            – index of the parameter that represents the initial time. 0 if fixed at 0.0.\n\nA mutable struct which represents a DSSCT model (see Discuit.jl models for further details).\n\n\n\n\n\n","category":"type"},{"location":"manual/#DPOMPs.SimResults","page":"DPOMPs.jl manual","title":"DPOMPs.SimResults","text":"SimResults\n\nFields\n\nmodel_name      – string, e,g, \"SIR\".\nparticle        – the 'trajectory' variable, of type Particle.\npopulation      – records the final system state.\nobservations    – simulated observations data (an Array of Observation types.)\n\nThe results of a simulation.\n\n\n\n\n\n","category":"type"},{"location":"manual/#Functions","page":"DPOMPs.jl manual","title":"Functions","text":"","category":"section"},{"location":"manual/#Simulation","page":"DPOMPs.jl manual","title":"Simulation","text":"","category":"section"},{"location":"manual/","page":"DPOMPs.jl manual","title":"DPOMPs.jl manual","text":"gillespie_sim","category":"page"},{"location":"manual/#DPOMPs.gillespie_sim","page":"DPOMPs.jl manual","title":"DPOMPs.gillespie_sim","text":"gillespie_sim(model, parameters; tmax = 100.0, num_obs = 5)\n\nRun a Doob-Gillespie simulation on model. Returns a SimResults type containing the trajectory and observations data, or an array of the same if n_sims > 1.\n\nParameters\n\nmodel       – DPOMPModel (see [DCTMPs.jl models]@ref).\nparameters  – model parameters.\n\nOptional\n\ntmax        – maximum time (default: 100.)\nn_obs       – the number of observations to draw (default: 5.)\nn_sims      – number of simulations to draw (default: 1.)\n\n\n\n\n\n","category":"function"},{"location":"manual/#Inference","page":"DPOMPs.jl manual","title":"Inference","text":"","category":"section"},{"location":"manual/","page":"DPOMPs.jl manual","title":"DPOMPs.jl manual","text":"run_mcmc_analysis\nrun_mbp_ibis_analysis\nrun_smc2_analysis\nrun_arq_mcmc_analysis","category":"page"},{"location":"manual/#DPOMPs.run_mcmc_analysis","page":"DPOMPs.jl manual","title":"DPOMPs.run_mcmc_analysis","text":"run_mcmc_analysis(model, obs_data, initial_parameters, steps = 50000, adapt_period = 10000, mbp = true, ppp = 0.3)\n\nParameters\n\nmodel               – DPOMPModel (see [DCTMPs.jl models]@ref).\nobs_data            – Observations data.\n\nOptional\n\nn_chains            – number of Markov chains (optional, default: 3.)\ninitial_parameters  – 2d array of initial model parameters. Each column vector correspondes to a single model parameter.\nsteps               – number of iterations.\nadapt_period        – number of discarded samples.\nmbp                 – model based proposals (MBP). Set mbp = false for standard proposals.\nppp                 – the proportion of parameter (vs. trajectory) proposals in Gibbs sampler. Default: 30%. NB. not required for MBP.\nfin_adapt           – finite adaptive algorithm. The default is false, i.e. [fully] adaptive.\n\nRun an n_chains-MCMC analysis. The initial_parameters are sampled from the prior distribution unless otherwise specified by the user.\n\nA Gelman-Rubin convergence diagnostic is automatically carried out for n_chains > 1 and included in the [multi-chain] results.\n\nOtherwise the results of a single-chain analysis are returned, which include the Geweke test statistics computed for that analysis.\n\n\n\n\n\n","category":"function"},{"location":"manual/#DPOMPs.run_mbp_ibis_analysis","page":"DPOMPs.jl manual","title":"DPOMPs.run_mbp_ibis_analysis","text":"run_mbp_ibis_analysis(model, obs_data, initial_parameters, ess_rs_crit = 0.5; n_props = 3, ind_prop = false, alpha = 1.002)\n\nRun an MBP IBIS analysis based on model and obs_data of type Observations.\n\nParameters\n\nmodel               – DPOMPModel (see [DCTMPs.jl models]@ref).\nobs_data            – Observations data.\nnp                  – number of particles (default = 2000.)\ness_rs_crit         – resampling criteria (default = 0.5.)\nn_props             – MBP mutations per step (default = 3.)\nind_prop            – true for independent theta proposals (default = false.)\nalpha               – user-defined, increase for lower acceptance rate targeted (default = 1.002.)\n\n\n\n\n\n","category":"function"},{"location":"manual/#DPOMPs.run_smc2_analysis","page":"DPOMPs.jl manual","title":"DPOMPs.run_smc2_analysis","text":"run_smc2_analysis(model, obs_data, initial_parameters, ess_rs_crit = 0.5; n_props = 3, ind_prop = false, alpha = 1.002)\n\nRun an SMC^2 (i.e. particle filter IBIS) analysis based on model and obs_data of type Observations.\n\nParameters\n\nmodel               – DPOMPModel (see [DCTMPs.jl models]@ref).\nobs_data            – Observations data.\nnp                  – number of particles (default = 2000.)\ness_rs_crit         – resampling criteria (default = 0.5.)\nn_props             – MBP mutations per step (default = 3.)\nind_prop            – true for independent theta proposals (default = false.)\nalpha               – user-defined, increase for lower acceptance rate targeted (default = 1.002.)\n\n\n\n\n\n","category":"function"},{"location":"manual/#Utilities","page":"DPOMPs.jl manual","title":"Utilities","text":"","category":"section"},{"location":"manual/#Visualisation","page":"DPOMPs.jl manual","title":"Visualisation","text":"","category":"section"},{"location":"manual/","page":"DPOMPs.jl manual","title":"DPOMPs.jl manual","text":"plot_trajectory\nplot_parameter_trace\nplot_parameter_marginal\nplot_parameter_heatmap","category":"page"},{"location":"manual/#DPOMPs.plot_trajectory","page":"DPOMPs.jl manual","title":"DPOMPs.plot_trajectory","text":"plot_trajectory(x)\n\nParameters\n\nx       – SimResults, i.e. from a call to gillespie_sim.\n\nPlot the trajectory of a a DGA simulation on model using UnicodePlots.jl.\n\n\n\n\n\n","category":"function"},{"location":"manual/#DPOMPs.plot_parameter_trace","page":"DPOMPs.jl manual","title":"DPOMPs.plot_parameter_trace","text":"plot_parameter_trace(mcmc, parameter)\n\nParameters\n\nsample      – MCMCSample, ARQMCMCSample or RejectionSample e.g. from a call to ADD XREF.\nparameter   – the index of the model parameter to be plotted.\n\nTrace plot of samples from n MCMC analyses for a given model parameter using UnicodePlots.jl.\n\n\n\n\n\n","category":"function"},{"location":"manual/#DPOMPs.plot_parameter_marginal","page":"DPOMPs.jl manual","title":"DPOMPs.plot_parameter_marginal","text":"plot_parameter_marginal(sample, parameter)\n\nPlot the marginal distribution of samples from an MCMC analysis for a given model parameter using UnicodePlots.jl.\n\nParameters\n\nresults     – Results object, e.g. of type MCMCSample.\nparameter   – the index of the model parameter to be plotted.\nadapt_period– Adaptation period to be discarded, only required for RejectionSample.\n\nOptional\n\nuse_is      – Resample IS rather than using MCMC [re]samples (ARQMCMCSample results only.)\n\n\n\n\n\n","category":"function"},{"location":"manual/#DPOMPs.plot_parameter_heatmap","page":"DPOMPs.jl manual","title":"DPOMPs.plot_parameter_heatmap","text":"plot_parameter_heatmap(mcmc, x_parameter, y_parameter)\n\nParameters\n\nmcmc        – MCMCResults, e.g. from a call to run_met_hastings_mcmc.\nx_parameter   – the index of the model parameter to be plotted on the x axis.\ny_parameter   – the index of the model parameter to be plotted on the y axis.\n\nPlot the marginal distribution of samples from an MCMC analysis for two model parameters using UnicodePlots.jl.\n\n\n\n\n\n","category":"function"},{"location":"manual/#Custom","page":"DPOMPs.jl manual","title":"Custom","text":"","category":"section"},{"location":"manual/#Index","page":"DPOMPs.jl manual","title":"Index","text":"","category":"section"},{"location":"manual/","page":"DPOMPs.jl manual","title":"DPOMPs.jl manual","text":"","category":"page"},{"location":"manual/#References","page":"DPOMPs.jl manual","title":"References","text":"","category":"section"},{"location":"manual/","page":"DPOMPs.jl manual","title":"DPOMPs.jl manual","text":"TBA","category":"page"},{"location":"#DPOMPs.jl","page":"DPOMPs.jl","title":"DPOMPs.jl","text":"","category":"section"},{"location":"","page":"DPOMPs.jl","title":"DPOMPs.jl","text":"Fast Bayesian parameter inference for Discrete-state-space Partially Observed Markov Processes in Julia","category":"page"}]
}
