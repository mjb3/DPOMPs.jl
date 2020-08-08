
## trajectory
"""
    plot_trajectory(x)

**Parameters**
- `x`       -- `SimResults`, i.e. from a call to `gillespie_sim`.

Plot the trajectory of a a DGA simulation on `model` using [UnicodePlots.jl](https://github.com/Evizero/UnicodePlots.jl).
"""
function plot_trajectory(x::SimResults)
    ## collect time and population
    t = zeros(length(x.particle.trajectory) + 1)
    # pop = x.population
    pop = zeros(Int64, length(x.particle.trajectory) + 1, length(x.particle.initial_condition))
    pop[1,:] .= x.particle.initial_condition
    for i in eachindex(x.particle.trajectory)
        t[i+1] = x.particle.trajectory[i].time
        pop[i+1, :] .= x.population[i]
    end
    ## plot
    p = UnicodePlots.lineplot(t, pop[:,1], title = string(x.model_name, " simulation"), name = string(x.model_name[1]))
    for i in 2:size(pop, 2)
        UnicodePlots.lineplot!(p, t, pop[:,i], name = string(x.model_name[i]))
    end
    UnicodePlots.xlabel!(p, "time")
    UnicodePlots.ylabel!(p, "population")
    return p
end

##
"""
    plot_parameter_trace(mcmc, parameter)

**Parameters**
- `sample`      -- `MCMCSample`, `ARQMCMCSample` or `RejectionSample` e.g. from a call to `ADD XREF`.
- `parameter`   -- the index of the model parameter to be plotted.

Trace plot of samples from `n` MCMC analyses for a given model `parameter` using [UnicodePlots.jl](https://github.com/Evizero/UnicodePlots.jl).
"""
function plot_parameter_trace(sample::RejectionSample, parameter::Int64)
    x = 1:size(sample.theta, 2)
    p = UnicodePlots.lineplot(x, sample.theta[parameter,:,1], title = string("θ", Char(8320 + parameter), " traceplot."))
    for i in 2:size(sample.theta, 3)
        UnicodePlots.lineplot!(p, sample.theta[parameter,:,i])
    end
    UnicodePlots.xlabel!(p, "sample")
    UnicodePlots.ylabel!(p, string("θ", Char(8320 + parameter)))
    return p
end

function plot_parameter_trace(sample::ARQMCMCSample, parameter::Int64)
    return plot_parameter_trace(sample.samples, parameter)
end

function plot_parameter_trace(sample::MCMCSample, parameter::Int64)
    return plot_parameter_trace(sample.samples, parameter)
end

## marginal
"""
    plot_parameter_marginal(sample, parameter)

Plot the marginal distribution of samples from an MCMC analysis for a given model `parameter` using [UnicodePlots.jl](https://github.com/Evizero/UnicodePlots.jl).

**Parameters**
- `results`     -- Results object, e.g. of type `MCMCSample`.
- `parameter`   -- the index of the model parameter to be plotted.
- `adapt_period`-- Adaptation period to be discarded, only required for `RejectionSample`.
**Optional**
- `use_is`      -- Resample IS rather than using MCMC [re]samples (`ARQMCMCSample` results only.)

"""
function plot_parameter_marginal(sample::RejectionSample, parameter::Int64, adapt_period::Int64)
    x = sample.theta[parameter, (adapt_period+1):size(sample.theta, 2), :]
    x = reshape(x, length(x))
    p = UnicodePlots.histogram(x, nbins = 20)
    UnicodePlots.ylabel!(p, string("θ", Char(8320 + parameter)))
    UnicodePlots.xlabel!(p, "samples")
    return p
end

## MCMC
function plot_parameter_marginal(sample::MCMCSample, parameter::Int64)
    return plot_parameter_marginal(sample.samples, parameter, sample.adapt_period)
end

## resampler
function plot_parameter_marginal(sample::ImportanceSample, parameter::Int64)
    rs = resample_is(sample)
    return plot_parameter_marginal(rs, parameter, 0)
end

function plot_parameter_marginal(sample::ARQMCMCSample, parameter::Int64; use_is::Bool = false)
    use_is && (return plot_parameter_marginal(sample.imp_sample, parameter))
    return plot_parameter_marginal(sample.samples, parameter, sample.adapt_period)
end


## heatmap
"""
    plot_parameter_heatmap(mcmc, x_parameter, y_parameter)

**Parameters**
- `mcmc`        -- `MCMCResults`, e.g. from a call to `run_met_hastings_mcmc`.
- `x_parameter`   -- the index of the model parameter to be plotted on the x axis.
- `y_parameter`   -- the index of the model parameter to be plotted on the y axis.

Plot the marginal distribution of samples from an MCMC analysis for two model parameters using [UnicodePlots.jl](https://github.com/Evizero/UnicodePlots.jl).
"""
function plot_parameter_heatmap(sample::RejectionSample, x_parameter::Int64, y_parameter::Int64, adapt_period::Int64)
    x = vec(sample.theta[x_parameter, (adapt_period+1):size(sample.theta,2), :])
    y = vec(sample.theta[y_parameter, (adapt_period+1):size(sample.theta,2), :])
    p = UnicodePlots.densityplot(x, y, color = :red)
    UnicodePlots.xlabel!(p, string("θ", Char(8320 + x_parameter)))
    UnicodePlots.ylabel!(p, string("θ", Char(8320 + y_parameter)))
    return p
end

## MCMC
function plot_parameter_heatmap(sample::MCMCSample, x_parameter::Int64, y_parameter::Int64)
    return plot_parameter_heatmap(sample.samples, x_parameter, y_parameter, sample.adapt_period)
end

## resampler
function plot_parameter_heatmap(sample::ImportanceSample, x_parameter::Int64, y_parameter::Int64)
    rs = resample_is(sample)
    return plot_parameter_heatmap(rs, x_parameter, y_parameter, 0)
end

## ARQ
function plot_parameter_heatmap(sample::ARQMCMCSample, x_parameter::Int64, y_parameter::Int64; use_is::Bool = false)
    use_is && (return plot_parameter_heatmap(sample.imp_sample, x_parameter, y_parameter))
    return plot_parameter_heatmap(sample.samples, x_parameter, y_parameter, sample.adapt_period)
end
