### for internal use:

### compute importance sample mean

## collect theta and weights
function collect_theta_weight(grid::Dict{Any,Any}, np::Int64)
    theta = zeros(np, length(grid))
    w = zeros(length(grid))
    for (index, value) in enumerate(grid)
        theta[:,index] .= value.second.sample
        w[index] = exp(value.second.log_likelihood)
    end
    return (theta, w) # ADD LENGTH OF GRID ***********
end

## compute autocorrelation for a single lag (C samples given mu and wcv)
function compute_autocorrelation(samples::Array{Float64, 2}, mu_var::Array{Float64, 2}, lag::Int64)
    output = zeros(size(mu_var, 1))
    # for each param:
    for j in eachindex(output)
        for i in 1:(size(samples, 1) - lag)
            output[j] += (samples[i,j] - mu_var[j,1]) * (samples[i + lag, j] - mu_var[j,1])
        end
        output[j] /= (size(samples, 1) - lag) * mu_var[j,2]
    end
    return output
end

## see Pooley 2018
function estimate_model_evidence(p_y::Float64)
    return -2 * log(p_y)
end


## print autocorrelation
"""
    print_autocorrelation(autocorrelation, fpath)

**Parameters**
- `autocorrelation` -- the results of a call to `compute_autocorrelation`.
- `fpath`           -- the file path of the destination file.

Save the results from a call to `compute_autocorrelation` to the file `fpath`, e.g. "./out/ac.csv".
"""
function print_autocorrelation(acr::AutocorrelationResults, fpath::String)
    open(fpath, "w") do f
        # print headers
        write(f, "lag")
        for j in 1:size(acr.autocorrelation, 2)
            write(f, ", x$j")
        end
        # print autocorr
        for i in 1:size(acr.autocorrelation, 1)
            # write(f, "\n$((i - 1) * AC_LAG_INT)")
            write(f, "\n$(acr.lag[i])")
            for j in 1:size(acr.autocorrelation, 2)
                write(f, ",$(acr.autocorrelation[i,j])")
            end
        end
    end
end
