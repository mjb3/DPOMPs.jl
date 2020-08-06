
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
