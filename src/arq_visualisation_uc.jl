
## traceplot
# """
#     plot_parameter_trace(mcmc, parameter)
#
# **Parameters**
# - `results`     -- `ARQMCMCResults`, i.e. from a call to `run_arq_mcmc_analysis`.
# - `parameter`   -- the index of the model parameter to be plotted.
#
# Trace plot of samples from `n` MCMC analyses for a given model `parameter` using [UnicodePlots.jl](https://github.com/Evizero/UnicodePlots.jl).
# """
# function plot_parameter_trace(results::ARQMCMCResults, parameter::Int64)
#     mcmc = results.mcmc
#     x = 1:size(mcmc[1].samples, 1)
#     p = UnicodePlots.lineplot(x, mcmc[1].samples[:, parameter], title = string("θ", Char(8320 + parameter), " traceplot."))
#     for i in 2:length(mcmc)
#         UnicodePlots.lineplot!(p, mcmc[i].samples[:, parameter])
#     end
#     UnicodePlots.xlabel!(p, "sample")
#     UnicodePlots.ylabel!(p, string("θ", Char(8320 + parameter)))
#     return p
# end

## marginal
# """
#     plot_parameter_marginal(mcmc, parameter)
#
# **Parameters**
# - `results`     -- `ARQMCMCResults`, i.e. from a call to `run_arq_mcmc_analysis`.
# - `parameter`   -- the index of the model parameter to be plotted.
#
# Plot the marginal distribution of samples from an MCMC analysis for a given model `parameter` using [UnicodePlots.jl](https://github.com/Evizero/UnicodePlots.jl).
# """
# function plot_parameter_marginal(results::ARQMCMCResults, parameter::Int64)
#     x = zeros(length(results.mcmc) * (size(results.mcmc[1].samples, 1) - results.adapt_period))
#     for i in 1:length(results.mcmc)
#         st = ((i - 1) * (size(results.mcmc[i].samples, 1) - results.adapt_period)) + 1
#         fn = i * (size(results.mcmc[i].samples, 1) - results.adapt_period)
#         x[st:fn] .= results.mcmc[i].samples[(results.adapt_period + 1):size(results.mcmc[i].samples, 1), parameter]
#     end
#     # x = Float64.(x)
#     # println("x type:", typeof(x))
#     p = UnicodePlots.histogram(x, nbins = results.grid_resolution)
#     UnicodePlots.ylabel!(p, string("θ", Char(8320 + parameter)))
#     UnicodePlots.xlabel!(p, "samples")
#     return p
# end

## heatmap
# """
#     plot_parameter_heatmap(mcmc, x_parameter, y_parameter)
#
# **Parameters**
# - `mcmc`        -- `MCMCResults`, e.g. from a call to `run_met_hastings_mcmc`.
# - `x_parameter`   -- the index of the model parameter to be plotted on the x axis.
# - `y_parameter`   -- the index of the model parameter to be plotted on the y axis.
#
# Plot the marginal distribution of samples from an MCMC analysis for two model parameters using [UnicodePlots.jl](https://github.com/Evizero/UnicodePlots.jl).
# """
# function plot_parameter_heatmap(results::ARQMCMCResults, x_parameter::Int64, y_parameter::Int64)
#     x = zeros(length(results.mcmc) * (size(results.mcmc[1].samples, 1) - results.adapt_period))
#     y = zeros(length(results.mcmc) * (size(results.mcmc[1].samples, 1) - results.adapt_period))
#     for i in 1:length(results.mcmc)
#         st = ((i - 1) * (size(results.mcmc[i].samples, 1) - results.adapt_period)) + 1
#         fn = i * (size(results.mcmc[i].samples, 1) - results.adapt_period)
#         x[st:fn] .= results.mcmc[i].samples[(results.adapt_period + 1):size(results.mcmc[i].samples, 1), x_parameter]
#         y[st:fn] .= results.mcmc[i].samples[(results.adapt_period + 1):size(results.mcmc[i].samples, 1), y_parameter]
#     end
#     p = UnicodePlots.densityplot(x, y, color = :red)
#     UnicodePlots.xlabel!(p, string("θ", Char(8320 + x_parameter)))
#     UnicodePlots.ylabel!(p, string("θ", Char(8320 + y_parameter)))
#     return p
# end

## autocorrelation R
"""
    plot_autocorrelation(autocorrelation)

**Parameters**
- `autocorrelation`     -- The results from a call to `compute_autocorrelation`.

Plot autocorrelation for an MCMC analysis.
"""
function plot_autocorrelation(autocorrelation::AutocorrelationResults)
    # build y
    for i in eachindex(autocorrelation.autocorrelation)
        autocorrelation.autocorrelation[i] = max(autocorrelation.autocorrelation[i], 0)
    end
    # plot
    p = UnicodePlots.lineplot(autocorrelation.lag, autocorrelation.autocorrelation[:,1], title = string("θ autocorrelation"))
    for i in 2:size(autocorrelation.autocorrelation, 2)
        UnicodePlots.lineplot!(p, autocorrelation.lag, autocorrelation.autocorrelation[:,i])
    end
    UnicodePlots.xlabel!(p, "lag")
    UnicodePlots.ylabel!(p, "R")
    return p
end
