### for internal use:

## does what it says on the tin SOON TO BE OBSOLETE ****
# function gelman_diagnostic(mcmc::Array{MCMCResults,1}, theta_size::Int64, num_iter::Int64)
#     ## compute W; B; V
#     # collect means and variances
#     mce = Array{Float64, 2}(undef, length(mcmc), theta_size)
#     mcv = Array{Float64, 2}(undef, length(mcmc), theta_size)
#     for i in eachindex(mcmc)
#         mce[i,:] .= mcmc[i].mean
#         for j in 1:theta_size
#             mcv[i,j] = mcmc[i].covar[j,j]
#         end
#     end
#     # compute W, B
#     b = Array{Float64, 1}(undef, theta_size)
#     w = Array{Float64, 1}(undef, theta_size)
#     mu = Array{Float64, 1}(undef, theta_size)
#     co = Array{Float64, 1}(undef, theta_size)
#     v = Array{Float64, 1}(undef, theta_size)
#     for j in 1:theta_size
#         b[j] = num_iter * Statistics.cov(mce[:,j])
#         w[j] = Statistics.mean(mcv[:,j])
#         # mean of means and var of vars (used later)
#         mu[j] = Statistics.mean(mce[:,j])
#         co[j] = Statistics.cov(mcv[:,j])
#         # compute pooled variance
#         v[j] = w[j] * ((num_iter - 1) / num_iter) + b[j] * ((theta_size + 1) / (theta_size * num_iter))
#     end
#     vv_w = Array{Float64, 1}(undef, theta_size)   # var of vars (theta_ex, i.e. W)
#     vv_b = Array{Float64, 1}(undef, theta_size)   # var of vars (B)
#     mce2 = copy(mce)
#     mce2 .*= mce                                           # ev(theta)^2
#     cv_wb = Array{Float64, 1}(undef, theta_size)   # wb covar
#     for j in 1:theta_size
#         vv_w[j] = co[j] / length(mcmc)
#         vv_b[j] = (2 * b[j] * b[j]) / (length(mcmc) - 1)
#         cv_wb[j] = (num_iter / length(mcmc)) * (Statistics.cov(mcv[:,j], mce2[:,j]) - (2 * mu[j] * Statistics.cov(mcv[:,j], mce[:,j])))
#     end
#     # compute d; d_adj (var.V)
#     d = Array{Float64, 1}(undef, theta_size)
#     dd = Array{Float64, 1}(undef, theta_size)
#     atmp = num_iter - 1
#     btmp = 1 + (1 / length(mcmc))
#     for j in 1:theta_size
#         tmp = ((vv_w[j] * atmp * atmp) + (vv_b[j] * btmp * btmp) + (cv_wb[j] * 2 * atmp * btmp)) / (num_iter * num_iter)
#         # println(" ", tmp)
#         d[j] = (2 * v[j] * v[j]) / tmp
#         dd[j] = (d[j] + 3) / (d[j] + 1)
#     end
#     # compute scale reduction estimate
#     sre = Array{Float64, 1}(undef, theta_size)
#     sre_ll = Array{Float64, 1}(undef, theta_size)
#     sre_ul = Array{Float64, 1}(undef, theta_size)
#     for j in 1:theta_size
#         rr = (1 + (1 / length(mcmc))) * (1 / num_iter)  * (b[j] / w[j])
#         sre[j] = sqrt(dd[j] * (((num_iter - 1) / num_iter) + rr))
#         # F dist(nu1, nu2)
#         fdst = Distributions.FDist(length(mcmc) - 1, 2 * w[j] * w[j] / vv_w[j])
#         sre_ll[j] = sqrt(dd[j] * (((num_iter - 1) / num_iter) + Statistics.quantile(fdst, 0.025) * rr))
#         sre_ul[j] = sqrt(dd[j] * (((num_iter - 1) / num_iter) + Statistics.quantile(fdst, 0.975) * rr))
#     end
#     # return results
#     return (mu, sqrt.(b), sqrt.(w), sre, sre_ll, sre_ul)
# end

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

## std grid, i.e. for each GridPoint TBR
# function compute_is_mean(grid::Dict{Any, Any}, np::Int64)
#     tpd = 0.0
#     mu = zeros(np)
#     for (key, val) in grid
#         pd = exp(val.log_likelihood)
#         tpd += pd
#         mu .+= val.sample * pd
#     end
#     mu ./= tpd
#     tpd /= length(grid)
#     return (is_mu, tpd) # ADD LENGTH OF GRID ***********
# end

## da grid, i.e. for each GridSample, within each GridSet
# function compute_da_is_mean(grid::Dict{Any, Any}, np::Int64)
#     tpd = 0.0
#     is_mu = zeros(np)
#     for (key, val) in grid
#         # val: ss
#         for s in eachindex(val.samples)
#             pd = exp(val.samples[s].log_likelihood)
#             tpd += pd
#             is_mu .+= (val.samples[s].sample * pd)
#         end
#     end
#     is_mu ./= tpd
#     return (is_mu, tpd) # ADD LENGTH OF GRID ***********
# end

### public utils

## autocorrelation R'
# """
#     compute_autocorrelation(mcmc, lags = 200)
#
# **Parameters**
# - `mcmc`    -- an array of `MCMCResults` variables.
# - `lags`    -- the number of lags to compute. Default: 200.
#
# Compute autocorrelation R' for an ARQMCMC analysis. The formula for multiple Markov chains is given by:
#
# \$R^{\\prime}_l = \\frac{\\textrm{E} [ (X_i - \\bar{X}_b) ( X_{i + l} - \\bar{X}_b ) ]}{\\sigma^2_b}\$
#
# \$\\sigma^2_b = \\textrm{E} [(X_i - \\bar{X}_b)^2]\$
#
# for any given lag `l` up to `lags` (default: 200).
# """
# function compute_autocorrelation(results::ARQMCMCResults, lags::Int64 = 200)
#     AC_LAG_INT = 10
#     ## NEW
#     mu_wcv = compute_mean_wcv(results.mcmc, results.adapt_period)
#     # mu = mu_wcv[:,1]
#     # wcv = mu_wcv[:,2]
#     ## for each lag interval
#     lag = zeros(Int64, lags + 1)
#     output = zeros(lags + 1, length(results.mu))
#     ## NEW
#     output[1,:] .= 1
#     for l in 1:lags
#         lag[l + 1] = l * AC_LAG_INT
#         for mc in eachindex(results.mcmc)
#             adp_s = results.mcmc[mc].samples[(results.adapt_period + 1):size(results.mcmc[mc].samples, 1), :]
#             output[l + 1,:] .+= compute_autocorrelation(adp_s, mu_wcv, lag[l + 1])
#         end
#         output[l + 1,:] ./= length(results.mcmc)
#     end
#     ## OLD (slightly more efficient but messy)
#     # for l in 0:lags
#     #     ## REINDEX HERE? ***********
#     #     lag[l + 1] = l * AC_LAG_INT
#     #     for j in eachindex(mu)
#     #         for mc in eachindex(mcmc)
#     #             for i in (results.adapt_period + 1):(size(mcmc[mc].samples, 1) - lag[l + 1])
#     #                 output[l + 1,j] += (mcmc[mc].samples[i,j] - mu[j]) * (mcmc[mc].samples[i + lag[l + 1], j] - mu[j])
#     #             end
#     #         end
#     #         output[l + 1,j] /= length(mcmc) * (size(mcmc[1].samples, 1) - results.adapt_period)
#     #         output[l + 1,j] /= wcv[j]
#     #     end
#     # end
#     return AutocorrelationResults(lag, output)
# end

## compute mean and 'whole chain' variance for n (adapted) chains
# function compute_mean_wcv(mcmc::Array{MCMCResults, 1}, adapt_period::Int64)
#     output = zeros(length(mcmc[1].mean), 2)
#     # collect means
#     mce = Array{Float64, 2}(undef, length(mcmc), length(mcmc[1].mean))
#     for mc in eachindex(mcmc)
#         mce[mc,:] .= mcmc[mc].mean
#     end
#     # for each parameter
#     for j in 1:size(output, 1)
#         # compute mean and whole chain var for adapted chains
#         output[j,1] = Statistics.mean(mce[:,j])
#         for mc in eachindex(mcmc)
#             for i in (adapt_period + 1):(size(mcmc[mc].samples, 1))
#                 output[j,2] += (mcmc[mc].samples[i,j] - output[j,1]) * (mcmc[mc].samples[i,j] - output[j,1])
#             end
#         end
#         output[j,2] /= length(mcmc) * (size(mcmc[1].samples, 1) - adapt_period)
#     end
#     return output
# end

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

## print individual MCMC results (internal use only)
# function print_arq_mcmc_result(mcmc::MCMCResults, dpath::String)
#     # create directory
#     isdir(dpath) || mkpath(dpath)
#     # print metadata (CONVERT TO ONE FILE FOR ALL CHAINS *)
#     open(string(dpath, "metadata.csv"), "w") do f
#         # print headers
#         write(f, "process_count")
#         # print md
#         write(f, "\n$(mcmc.process_count))")
#     end
#     # print parameter summary
#     open(string(dpath, "parameters.csv"), "w") do f
#         # print headers
#         write(f, "parameter,mean,sd")
#         for p in eachindex(mcmc.mean)
#             sd = sqrt(mcmc.covar[p,p])
#             write(f, "\n$p,$(mcmc.mean[p]),$sd")
#         end
#     end
#     ## print samples
#     open(string(dpath, "samples.csv"), "w") do f
#         # print headers
#         write(f, "iter")
#         for p in 1:size(mcmc.samples, 2)
#             write(f, ",x$p")
#         end
#         #
#         for i in 1:size(mcmc.samples, 1)
#             write(f, "\n$i")
#             for p in 1:size(mcmc.samples, 2)
#                 # t = mc[i, p]
#                 write(f, ",$(mcmc.samples[i, p])")
#             end
#         end
#     end # end of print samples
#     ## print sample metadata
#     open(string(dpath, "sample_md.csv"), "w") do f
#         # print headers
#         write(f, "iter,accepted,xf_ll,sys_time")
#         for p in 1:size(mcmc.samples, 2)
#             write(f, ",xi$p,xf$p")
#         end
#         # print metadata
#         for i in 1:size(mcmc.samples, 1)
#             write(f, "\n$i,$(mcmc.mc_accepted[i]),$(mcmc.mc_log_like[i]),$(mcmc.mc_time[i])")
#             for p in 1:size(mcmc.samples, 2)
#                 # t = mc[i, p]
#                 write(f, ",$(mcmc.sample_idx[i, p]),$(mcmc.mcf[i, p])")
#             end
#         end
#     end # end of print samples
# end
#
# ## print std grid
# function print_imp_sample(results::ARQMCMCResults, fpath::String)
#     open(fpath, "w") do f
#         # print headers
#         write(f, "visited,sampled,log_like")
#         for p in eachindex(results.mu)
#             write(f, ", i$p, p$p")
#         end
#         # print grid
#         for (key, val) in results.grid
#             write(f, "\n$(val.visited),$(val.sampled),$(val.log_likelihood)")
#             for p in eachindex(results.mu)
#                 write(f, ",$(key[p]),$(val.sample[p])")
#             end
#        end
#    end  # end of print grid
# end
#
# ## print DA grid
# # i.e. results.algorithm == C_ALG_STD?
# function print_da_imp_sample(results::ARQMCMCResults, fpath::String)
#     open(fpath, "w") do f
#         # print headers
#         write(f, "i,log_like")
#         for p in eachindex(results.mu)
#             write(f, ", i$p, p$p")
#         end
#         # print grid
#         i = 0
#         for (key, val) in results.grid
#             i += 1
#             for s in eachindex(val.samples)
#                 write(f, "\n$i,$(val.samples[s].log_likelihood)")
#                 for p in eachindex(results.mu)
#                     write(f, ",$(key[p]),$(val.samples[s].sample[p])")
#                 end
#             end
#        end
#    end  # end of print grid
# end

# ## DA grid sample
# struct GridSample
#     sample::Array{Float64, 1}   # i.e. theta
#     log_likelihood::Float64
# end
#
# ## DA grid 'value'
# struct GridSet
#     anchor::Array{Float64, 1}   # i.e. theta
#     samples::Array{GridSample, 1}
#     # log_likelihood::Float64     # i.e weighted density
# end

## PRINT ARQMCMC RESULTS
# """
#     print_arq_mcmc_results
#
# **Fields**
# - `results`     -- `ARQMCMCResults` (i.e. the output from an analysis).
# - `dpath`       -- the path of the directory where the output will be saved.
#
# Print the results of a limited sample size quasi MCMC (ARQMCMC) analysis in the desired directory.
# """
# function print_arq_mcmc_results(results::ARQMCMCResults, dpath::String)
#     print("printing: ", dpath)
#     # create directory
#     # dpath = string("./out/", dname, "/")
#     isdir(dpath) || mkpath(dpath)
#     # print metadata NEED TO FINISH ************
#     open(string(dpath, "metadata.csv"), "w") do f
#         # print headers
#         write(f, "algorithm,grid_res,chains,params,sample_limit,jitter,iterations,adapt_period,is_nc") # process_count, run_time, prop_type
#         nmc = length(results.mcmc)
#         # print md
#         write(f, "\n$(results.algorithm),$(results.grid_resolution),$nmc,$(length(results.mu)),$(results.sample_limit),$(results.jitter),$(results.iterations),$(results.adapt_period),$(results.is_bme)")
#     end
#     # print summary by theta row
#     open(string(dpath, "parameters.csv"), "w") do f
#         # print headers
#         write(f, "parameter,mu,is_mu,sdb,sdw,rng_ll,rng_ul,sre,sre_ll,sre_ul")
#         for p in eachindex(results.mu)
#             write(f, "\n$p,$(results.mu[p]),$(results.is_mu[p]),$(results.sdb[p]),$(results.sdw[p]),$(results.grid_range[p,1]),$(results.grid_range[p,2]),$(results.sre[p]),$(results.sre_ll[p]),$(results.sre_ul[p])")
#         end
#     end # end of print summary
#     # print chains
#     for i in eachindex(results.mcmc)
#         print_arq_mcmc_result(results.mcmc[i], string(dpath, "mc", i, "/"))
#     end
#     ## PRINT GRID
#     if results.algorithm == C_ALG_STD
#         print_imp_sample(results, string(dpath, "grid.csv"))
#     else
#         print_da_imp_sample(results, string(dpath, "grid.csv"))
#     end
#     println(" - done.")
# end

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

### tabulate stuff

# ## proposal summary
# function tabulate_proposals(results::ARQMCMCResults)
#     println("Proposal summary ", length(results.mcmc)," chains):")
#     h = ["Adapted", "Proposed", "Accepted", "Rate", "f(x)"]
#     ## summary
#     pr = [length(results.mcmc) * results.adapt_period, length(results.mcmc) * (results.iterations - results.adapt_period)]
#     ac = zeros(Int64, 2)
#     fx = zeros(Int64, 2)
#     for mc in eachindex(results.mcmc)
#         ac[1] += sum(results.mcmc[mc].mc_accepted[1:results.adapt_period])
#         ac[2] += sum(results.mcmc[mc].mc_accepted[(results.adapt_period + 1):results.iterations])
#         fx[1] += sum(results.mcmc[mc].mc_fx[1:results.adapt_period])
#         fx[2] += sum(results.mcmc[mc].mc_fx[(results.adapt_period + 1):results.iterations])
#     end
#     d = Matrix(undef, 3, 5)
#     d[1, :] .= [ "false", pr[1], ac[1], round(100 * ac[1] / pr[1]; sigdigits = C_PR_SIGDIG), fx[1] ]
#     d[2, :] .= [ "true", pr[2], ac[2], round(100 * ac[2] / pr[2]; sigdigits = C_PR_SIGDIG), fx[2] ]
#     d[3, :] .= [ "total", sum(pr), sum(ac), round(100 * sum(ac) / sum(pr); sigdigits = C_PR_SIGDIG), sum(fx) ]
#     ## detailed TBR
#     # d = Matrix(undef, length(results.mcmc) * 2, 5)
#     # for mc in eachindex(results.mcmc)
#     #     # burnin
#     #     pr = results.adapt_period
#     #     ac = sum(results.mcmc[mc].mc_accepted[1:results.adapt_period])
#     #     d[mc, :] .= [ "false", mc, pr, ac, round(100 * ac / pr; sigdigits = C_PR_SIGDIG) ]
#     # end
#     # for mc in eachindex(results.mcmc)
#     #     # adapted
#     #     pr = length(results.mcmc[mc].mc_accepted) - results.adapt_period
#     #     ac = sum(results.mcmc[mc].mc_accepted[(results.adapt_period + 1):length(results.mcmc[mc].mc_accepted)])
#     #     d[(length(results.mcmc) + mc), :] .= [ "true", mc, pr, ac, round(100 * ac / pr; sigdigits = C_PR_SIGDIG) ]
#     # end
#     PrettyTables.pretty_table(d, h)
# end
#
# ## results summary AS DataFrame?
# """
#     tabulate_samples
#
# **Parameters**
# - `results`     -- the results of a call to `run_arq_mcmc_analysis`.
# - `proposals`   -- display proposal analysis.
#
# Display the results of an `ARQMCMC` analysis.
# """
# function tabulate_samples(results::ARQMCMCResults, proposals = false)
#     ## proposals
#     proposals && tabulate_proposals(results)
#     ## samples
#     println(results.algorithm, " MCMC results:")
#     h = ["θ", "Iμ", "Rμ", "σ", "SRE", "SRE95", "BME"]
#     d = Matrix(undef, length(results.mu), 7)
#     d[:,1] .= 1:length(results.mu)
#     d[:,2] .= round.(results.is_mu; sigdigits = C_PR_SIGDIG)
#     d[:,3] .= round.(results.mu; sigdigits = C_PR_SIGDIG)
#     d[:,4] .= round.(results.sdw; sigdigits = C_PR_SIGDIG)
#     d[:,5] .= round.(results.sre; sigdigits = C_PR_SIGDIG)
#     d[:,6] .= round.(results.sre_ul; sigdigits = C_PR_SIGDIG)
#     d[:,7] .= 0.0
#     d[1,7] = round(results.is_bme; sigdigits = C_PR_SIGDIG)
#     PrettyTables.pretty_table(d, h)
# end
