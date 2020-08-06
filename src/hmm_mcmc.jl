#### data augmented MCMC ####

## global const
const C_INITIAL = 0.1                    # proposal scalar

## proposal step
macro met_hastings_init()
    esc(quote
    steps::Int64 = size(theta, 2) #GET RID
    ADAPT_INTERVAL = adapt_period / 10   # interval between adaptation steps
    xi = x0
    ## covar matrix
    covar = zeros(length(xi.theta), length(xi.theta))
    for i in eachindex(xi.theta)
        covar[i,i] = xi.theta[i] == 0.0 ? 1 : (xi.theta[i]^2)
    end
    propd = Distributions.MvNormal(covar)
    c = C_INITIAL
    theta[:,1,mc] .= xi.theta               # samples
    end)
end

## adaptation
macro met_hastings_adapt()
    esc(quote
    ## adaptation
    if (!fin_adapt || i < adapt_period)
        c *= (accepted ? 1.002 : 0.999)
        # end of adaption period
        if i % ADAPT_INTERVAL == 0
            covar = Distributions.cov(transpose(theta[:,1:i,mc]))
            if sum(covar) == 0
                println("warning: low acceptance rate detected in adaptation period")
            else
                # C_DEBUG && println(" cv := ", transpose(theta[:,1:i,mc]))
                propd = Distributions.MvNormal(covar)
            end
        end
    end
    end)
end

## TO DO
# - break with macros
# - add Gibbs sampler
# - add std proposal example
# - add custom example for Robs

## Single particle adaptive Metropolis Hastings algorithm
function met_hastings_alg!(theta::Array{Float64,3}, mc::Int64, model::HiddenMarkovModel, adapt_period::Int64, x0::Particle, proposal_alg::Function, fin_adapt::Bool) #, joint_proposals::Bool, ppp::Float64
    @met_hastings_init
    for i in 2:steps            # met_hastings_step
        xf = proposal_alg(model, get_mv_param(propd, c, theta[:,i-1,mc]), xi)
        if sum(xf.log_like) == -Inf
            accepted = false    # reject automatically
        else                    # accept or reject
            # NB: [2] == full g(x) log like
            mh_prob::Float64 = exp(xf.log_like[2] - xi.log_like[2])   # HACK:
            accepted = (mh_prob > 1 || mh_prob > rand())
        end
        accepted && (xi = xf)
        theta[:,i,mc] .= xi.theta
        ## adaptation
        @met_hastings_adapt
    end ## end of MCMC loop
end

## Single particle adaptive Gibbs sampler - TO BE FINISHED ****
function gibbs_mh_alg!(theta::Array{Float64,3}, mc::Int64, model::HiddenMarkovModel, adapt_period::Int64, x0::Particle, proposal_alg::Function, fin_adapt::Bool, ppp::Float64)
    @met_hastings_init
    for i in 2:steps            # Gibbs
        xf = proposal_alg(model, get_mv_param(propd, c, theta[:,i-1,mc]), xi)
        if sum(xf.log_like) == -Inf
            accepted = false    # reject automatically
        else                    # accept or reject
            # NB: [2] == full g(x) log like
            mh_prob::Float64 = exp(xf.log_like[2] - xi.log_like[2])
            accepted = (mh_prob > 1 || mh_prob > rand())
        end
        accepted && (xi = xf)
        theta[:,i,mc] .= xi.theta
        ## adaptation
        @met_hastings_adapt
    end ## end of MCMC loop
end

## generic met hastings
# function generic_mcmc!(theta::Array{Float64,3}, mc::Int64, model::HiddenMarkovModel, adapt_period::Int64, theta0::Array{Float64,1}, target_log_density::Function) #, joint_proposals::Bool, ppp::Float64
#     steps::Int64 = size(theta,2) #GET RID
#     ADAPT_INTERVAL = adapt_period / 10   # interval between adaptation steps
#     ll_i = target_log_density(theta0)
#     # covar matrix
#     covar = zeros(length(theta0), length(theta0))
#     for i in eachindex(theta0)
#         covar[i,i] = 0.1 * (theta0[i] == 0.0 ? 1 : theta0[i]^2)
#     end
#     propd = Distributions.MvNormal(covar)
#     c = C_INITIAL
#     ## results
#     # theta = Array{Float64, 2}(undef, steps, length(theta0))
#     theta[:,mc,1,:] .= theta0 # FIX *
#     for i in 2:steps
#         theta[mc,i,:] = get_mv_param(propd, c, theta[mc,i-1,:])
#         ll_f = target_log_density(theta[mc,i,:])
#         if ll_f == -Inf
#             # reject automatically
#             accepted = false
#         else
#             # accept or reject
#             mh_prob::Float64 = exp(ll_f - ll_i)
#             accepted = (mh_prob > 1 || mh_prob > rand())
#         end
#         if accepted
#             ll_i = ll_f
#         else
#             theta[mc,i,:] .= theta[mc,i-1,:]
#         end
#         ## adaptation
#         if i < adapt_period
#             c *= (accepted ? 1.002 : 0.999)
#             # end of adaption period
#             if i % ADAPT_INTERVAL == 0
#                 covar = Distributions.cov(theta[:,1:i,mc])
#                 # output.cv .= Statistics.cov(reshape(theta[:,(ap+1):size(theta,3),:], d, size(theta,1)))
#                 if sum(covar) == 0
#                     println("warning: low acceptance rate detected in adaptation period")
#                 else
#                     propd = Distributions.MvNormal(covar)
#                 end
#             end
#         end
#     end ## end of MCMC loop
# end

## gelman diagnostic (internal)
function gelman_diagnostic(samples::Array{Float64,3}, discard::Int64)
    np = size(samples,1)
    niter::Int64  = size(samples,2)
    nmc = size(samples,3)
    fsmpl = discard + 1
    nsmpl = niter - discard
    ## compute W; B; V
    # collect means and variances
    mce = zeros(nmc, np)
    mcv = zeros(nmc, np)
    for i in 1:nmc
        for j in 1:np
            mce[i,j] = Statistics.mean(samples[j,fsmpl:end,i])
            mcv[i,j] = Statistics.cov(samples[j,fsmpl:end,i])
        end
    end
    # compute W, B
    b = zeros(np)
    w = zeros(np)
    mu = zeros(np)
    co = zeros(np)
    v = zeros(np)
    for j in 1:np
        b[j] = nsmpl * Statistics.cov(mce[:,j])
        w[j] = Statistics.mean(mcv[:,j])
        # mean of means and var of vars (used later)
        mu[j] = Statistics.mean(mce[:,j])
        co[j] = Statistics.cov(mcv[:,j])
        # compute pooled variance
        # - COMPARE THESE..?
        v[j] = w[j] * ((nsmpl - 1) / nsmpl) + b[j] * ((np + 1) / (np * nsmpl))
        # v[j] = ((num_s - 1) * (w[j] / num_s)) + ((1 + (1 / size(theta_init, 1))) * (b[j] / num_s))
        # v[j] = ((w[j] / num_s) * (num_s - 1)) + ((b[j] / num_s) * (1 + (1 / size(theta_init, 1))))
    end
    #
    vv_w = zeros(np)   # var of vars (theta_ex, i.e. W)
    vv_b = zeros(np)   # var of vars (B)
    mce2 = mce.^2                                          # ev(theta)^2
    cv_wb = zeros(np)   # wb covar
    for j in 1:np
        vv_w[j] = co[j] / nmc
        vv_b[j] = (2 * b[j] * b[j]) / (nmc - 1)
        cv_wb[j] = (nsmpl / nmc) * (Statistics.cov(mcv[:,j], mce2[:,j]) - (2 * mu[j] * Statistics.cov(mcv[:,j], mce[:,j])))
    end
    # compute d; d_adj (var.V)
    d = zeros(np)
    dd = zeros(np)
    atmp = nsmpl - 1
    btmp = 1 + (1 / nmc)
    for j in 1:np
        tmp = ((vv_w[j] * atmp * atmp) + (vv_b[j] * btmp * btmp) + (cv_wb[j] * 2 * atmp * btmp)) / (nsmpl * nsmpl)
        d[j] = (2 * v[j] * v[j]) / tmp
        dd[j] = (d[j] + 3) / (d[j] + 1)
    end
    # compute scale reduction estimate
    sre = zeros(np,3)
    try
        # sre_ll = zeros(np)
        # sre_ul = zeros(np)
        for j in 1:np
            rr = btmp * (1 / nsmpl)  * (b[j] / w[j]) ## NB. btmp ***
            sre[j,2] = sqrt(dd[j] * ((atmp / nsmpl) + rr)) ## atmp
            # F dist(nu1, nu2)
            fdst = Distributions.FDist(nmc - 1, 2 * w[j] * w[j] / vv_w[j])
            sre[j,1] = sqrt(dd[j] * ((atmp / nsmpl) + Statistics.quantile(fdst, 0.025) * rr))
            sre[j,3] = sqrt(dd[j] * ((atmp / nsmpl) + Statistics.quantile(fdst, 0.975) * rr))
        end
        # return GelmanResults
        return (mu = mu, wcv = sqrt.(w), sre = sre)
    catch gmn_err
        println("GELMAN ERROR: ", gmn_err)
        return (mu = mu, wcv = sqrt.(w), sre = sre) # return zeros
    end
end


#### public functions ####

## MBP MCMC 'analysis', i.e. gelman diagnostic
function run_mbp_mcmc(model::HiddenMarkovModel, theta_init::Array{Float64,2}, steps::Int64, adapt_period::Int64, fin_adapt::Bool)
    start_time = time_ns()
    # samples = zeros(size(theta_init,1), steps, size(theta_init,2))
    samples = zeros(size(theta_init,1), steps, size(theta_init,2))
    println("Running: ", size(theta_init, 2) ,"-chain ", steps, "-sample ", fin_adapt ? "finite-" : "", "adaptive MBP-MCMC analysis (model: ", model.model_name, ")")
    for i in 1:size(theta_init,2)
        x0 = generate_x0(model, theta_init[:,i])    # simulate initial particle
        ## run inference
        C_DEBUG && print( " mc", i, " initialising for x0 := ", x0.theta, " (", length(x0.trajectory), " events)")
        met_hastings_alg!(samples, i, model, adapt_period, x0, model_based_proposal, fin_adapt)
        println(" chain ", i, " complete.")
    end
    # rejs = RejectionSample(samples, zeros(size(theta_init,2)), zeros(size(theta_init,2), size(theta_init,2)))
    # compute_rs_mu_cv!(rejs.mu, rejs.cv, samples)
    rejs = handle_rej_samples(samples, adapt_period)
    gd = gelman_diagnostic(samples, adapt_period)         # run convergence diagnostic
    output = MCMCSample(rejs, adapt_period, gd.sre, time_ns() - start_time)
    # gelman_diagnostic!(output)                      # run convergence diagnostic
    println("- finished in ", Int64(round(output.run_time / C_RT_UNITS)), "s. E(x) := ", output.samples.mu)
    return output
end

## particle MCMC
# - NEED TO FIX THIS FOR NEW THETA LAYOUT
function run_pmcmc(model::HiddenMarkovModel, theta_init::Array{Float64,2}, steps::Int64 = 50000, adapt_period::Int64 = 10000, p::Int64 = 200)
    start_time = time_ns()
    samples = Array{Array{Float64,2},1}(undef, size(theta_init,2))
    println("Running PMCMC analysis: " , size(theta_init, 2), " x ", steps, " samples")
    ## target density
    ps = length(model.fn_initial_condition())
    function comp_log_pdf(theta::Array{Float64, 1})
        return model.fn_log_prior(theta) + estimate_likelihood(model, theta, p, ps, rsp_systematic)
    end
    ## run inference
    for i in 1:size(theta_init,1)
        generic_mcmc!(samples, i, model, adapt_period, theta_init[:,i], comp_log_pdf)
        println(" chain ", i, " complete.")
    end
    # rejs = RejectionSample(samples, zeros(size(theta_init,2)), zeros(size(theta_init,2)))
    # compute_rs_mu_cv!(rejs.mu, rejs.cv, samples)
    rejs = handle_rej_samples(samples, adapt_period)
    output = MCMCSample(rejs, adapt_period, zeros(size(theta_init,2), 3), time_ns() - start_time)
    gelman_diagnostic!(output)                      # run convergence diagnostic
    println("- finished in ", Int64(output.run_time / C_RT_UNITS), ". E(x) := ", output.samples.mu)
    return output
end
