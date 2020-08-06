## iterated batch importance sampling (in dev)

### resample algorithms (IBIS)
include("hmm_resample.jl")

## compute effective sample size
function compute_ess(w::Array{Float64,1})
    return sum(w)^2 / sum(w.^2)
end

## compute is mu var
function compute_is_mu_covar!(mu::Array{Float64,1}, cv::Array{Float64,2}, theta::Array{Float64,2}, w::Array{Float64,1})
    for i in eachindex(mu)
        mu[i] = sum(w .* theta[i,:]) / sum(w)                       # compute mu
        cv[i,i] = sum(w .* ((theta[i,:] .- mu[i]).^2)) / sum(w)     # variance
        for j in 1:(i-1)                                            # covar
            cv[i,j] = cv[j,i] = sum(w .* (theta[i,:] .- mu[i]) .* (theta[j,:] .- mu[j])) / sum(w)
        end
    end
end

## i.e. SMC^2 algorithm (Chopin 2013)
function run_pibis(model::HiddenMarkovModel, theta::Array{Float64, 2}, ess_rs_crit::Float64, ind_prop::Bool, alpha::Float64, np::Int64; n_props = 1)
    ##, rej_sample = size(theta, 2)
    outer_p = size(theta, 2)
    println("Running: ", outer_p, "-particle SMC^2 analysis (model: ", model.model_name, ")")
    start_time = time_ns()

    ess_crit = ess_rs_crit * outer_p
    fn_rs = rsp_systematic

    w = ones(outer_p)       # incremental weights g(y)
    aw = zeros(outer_p)     # ancestral weights (log likelihood)
    for i in eachindex(w)
        aw[i] = Distributions.logpdf(model.prior, theta[:,i])
    end
    ## initialise population matrices
    pop = Array{Array{Int64,2},1}(undef, outer_p)
    pop_size = length(model.fn_initial_condition())
    for i in 1:outer_p
        pop[i] = zeros(Int64, np, pop_size)
    end
    ## resampling workspace
    theta2 = copy(theta)
    aw2 = copy(aw)
    pop2 = deepcopy(pop)

    ## initialise and run:
    mu = zeros(size(theta, 1))
    cv = zeros(size(theta, 1), size(theta, 1))
    k_log = zeros(Int64, 2)
    bme = zeros(2)
    gx = zeros(outer_p)
    mtd_gx = zeros(outer_p)

    ## proposal distribution and scalar
    propd = Distributions.MvNormal(Matrix{Float64}(LinearAlgebra.I, size(theta,1), size(theta,1)))
    tj = 0.2

    ## for each observation
    obs_min = 1
    for obs_i in eachindex(model.obs_data)
        if model.obs_data[obs_i].obs_id > 0
            ## for each 'outer' particle
            for p in eachindex(pop)
                ## compute incremental weights (i.e. run pf)
                gx[p] = partial_log_likelihood!(pop[p], model, theta[:,p], fn_rs, ess_rs_crit, obs_min, obs_i)
            end
            aw .+= gx # update ancestral weights
            ## COMPUTE L
            gx .= exp.(gx)
            lml = log(sum(w .* gx) / sum(w))
            bme[1] += lml
            w .*= gx
            compute_is_mu_covar!(mu, cv, theta, w)
            ## resample and mutate if criteria satisfied:
            if compute_ess(w) < ess_crit
                C_DEBUG && println(" resampling for ", obs_i)
                ## update proposal density
                # propd = Distributions.MvNormal(cv)
                propd = get_prop_density(cv, propd)
                ## resample and swap
                nidx = rs_systematic(w)
                # pcnt = zeros(outer_p)   # abstract this ***
                for p in eachindex(pop)
                    # pcnt[nidx[p]] += 1  # vectorise this ***
                    theta2[:, p] .= theta[:, nidx[p]]
                    aw2[p] = aw[nidx[p]]
                    pop2[p] .= pop[nidx[p]]
                end
                mlr = Statistics.mean(gx[nidx]) * exp(lml)
                theta, theta2 = theta2, theta
                aw, aw2 = aw2, aw
                pop, pop2 = pop2, pop
                # - MUTATE:
                k_log[1] += outer_p
                mtd_gx .= gx[nidx]
                for p in eachindex(pop)
                    ## mutation step
                    for mki in 1:n_props
                        ## propose new theta - independent OR conditional on current sample
                        theta_f = ind_prop ? get_mv_param(propd, 1.0, mu) : get_mv_param(propd, tj, theta[:,p])
                        prtf = Distributions.logpdf(model.prior, theta_f)
                        if prtf != -Inf
                            pop_f = zeros(Int64, np, pop_size)
                            ## sample using pmcmc kernel
                            # - RETRIEVE MOST RECENT MARGINAL
                            if obs_i == 1
                                # - compute ml for obsi
                                gx_f = partial_log_likelihood!(pop_f, model, theta_f, fn_rs, ess_rs_crit, obs_i, obs_i)
                                aw_f = gx_f
                            else
                                # - compute aw up to obsi-1
                                aw_f = partial_log_likelihood!(pop_f, model, theta_f, fn_rs, ess_rs_crit, 1, obs_i - 1)
                                # - compute ml for obsi
                                gx_f = partial_log_likelihood!(pop_f, model, theta_f, fn_rs, ess_rs_crit, obs_i, obs_i)
                                aw_f += gx_f
                            end
                            # aw_f = prtf + partial_log_likelihood!(pop_f, model, theta_f, fn_rs, ess_rs_crit, 1, obs_i)
                            aw_f += prtf
                            ## accept with p(mh)
                            if exp(aw_f - aw[p]) > rand()
                                mtd_gx[p] = exp(gx_f)
                                theta[:,p] .= theta_f
                                aw[p] = aw_f
                                pop[p] .= pop_f
                                k_log[2] += 1
                                tj *= alpha
                            else
                                # mtd_gx[p] = mtd_gx[p]
                                tj *= 0.999
                            end
                        end
                    end
                end
                ## RB ML update
                bme[2] += log(mlr / Statistics.mean(mtd_gx))
                ## reset w = 1
                w .= 1
            else
                ## standard ML update
                bme[2] += log(sum(w .* gx) / sum(w))
            end # END OF sample/mutate
            obs_min = obs_i + 1
        end
    end # END OF OBS LOOP
    ## return weighted importance sample
    compute_is_mu_covar!(mu, cv, theta, w)
    C_DEBUG && println(" mcv: ", mu, cv)
    bme .*= -2
    output = ImportanceSample(mu, cv, theta, w, time_ns() - start_time, bme)
    println("- finished in ", Int64(round(output.run_time / C_RT_UNITS)), "s. AR = ", round(100.0 * k_log[2] / k_log[1]; sigdigits = 3), "%")
    return output
end

##
function get_prop_density(cv::Array{Float64,2}, old)
    ## update proposal density
    tmp = LinearAlgebra.Hermitian(cv)
    if LinearAlgebra.isposdef(tmp)
        return Distributions.MvNormal(Matrix(tmp))
    else
        println("warning: particle degeneracy problem")
        println(" covariance matrix: ", cv)
        return old
    end
end

## DEBUG:
# function copy_particle(p::Particle)
#     traj = Event[]
#     for i in eachindex(p.trajectory)
#         push!(traj, deepcopy(p.trajectory[i]))
#     end
#     return Particle(copy(p.theta), copy(p.initial_condition), copy(p.final_condition), traj, copy(p.log_like))
# end

## MBP IBIS algorithm
function run_mbp_ibis(model::HiddenMarkovModel, theta::Array{Float64, 2}, ess_rs_crit, n_props, ind_prop, alpha)
    outer_p = size(theta,2)
    println("Running: ", outer_p, "-particle MBP-IBIS analysis (model: ", model.model_name, ")")
    C_DEBUG && println(" - {n_props, ind_prop, ess} : = ", (n_props, ind_prop, ess_rs_crit))
    start_time = time_ns()
    ess_crit = ess_rs_crit * outer_p
    fn_rs = rsp_systematic

    ## initialise particles
    ptcls = Array{Particle,1}(undef, outer_p)
    # theta = copy(theta_init) # GET RID? *
    for p in eachindex(ptcls)
        ic = model.fn_initial_condition()
        ptcls[p] = Particle(theta[:,p], ic, copy(ic), Event[], [Distributions.logpdf(model.prior, theta[:,p]), 0.0])
    end
    ## resampling workspace
    ptcls2 = deepcopy(ptcls)

    ## proposal distribution and scalar
    propd = Distributions.MvNormal(Matrix{Float64}(LinearAlgebra.I, size(theta,1), size(theta,1)))
    tj = 0.2

    ## initialise and run:
    w = ones(outer_p)
    mu = zeros(size(theta, 1))
    cv = zeros(size(theta, 1), size(theta, 1))
    k_log = zeros(Int64, 2)
    ## for each observation
    bme = zeros(2)
    gx = zeros(outer_p)
    mtd_gx = zeros(outer_p)
    t = zeros(outer_p)
    model.t0_index > 0 && (t .= theta[model.t0_index,:])
    for obs_i in eachindex(model.obs_data)
        if model.obs_data[obs_i].obs_id > 0
            ## for each 'outer' particle
            for p in eachindex(ptcls)
                ## compute incremental weights (i.e. run pf)
                gx[p] = exp(iterate_particle!(ptcls[p], model, t[p], model.obs_data[obs_i]))
            end
            ## COMPUTE L and update weights
            lml = log(sum(w .* gx) / sum(w))
            bme[1] += lml
            # bme[1] += log(sum(w .* gx) / sum(w))
            w .*= gx
            ##
            compute_is_mu_covar!(mu, cv, theta, w)
            ## resample and mutate if criteria satisfied:
            essv = compute_ess(w)
            # t = model.obs_data[obs_i].time
            if (essv < ess_crit)
                C_DEBUG && println(" resampling for ", obs_i)
                ## resample and swap
                propd = get_prop_density(cv, propd)
                nidx = rs_systematic(w)
                # println(obs_i, " - ", length(nidx), " - ", length(mtd_gx))
                mtd_gx .= gx[nidx]
                for p in eachindex(ptcls)
                    ptcls2[p] = deepcopy(ptcls[nidx[p]])
                end
                mlr = Statistics.mean(gx[nidx]) * exp(lml)
                ptcls, ptcls2 = ptcls2, ptcls
                # mutate:
                k_log[1] += outer_p * n_props
                for p in eachindex(ptcls)
                    for mki in 1:n_props
                        ## propose new theta - independent OR conditional on current sample (rec'd)
                        theta_f = ind_prop ? get_mv_param(propd, 1.0, mu) : get_mv_param(propd, tj, ptcls[p].theta)
                        xf = partial_model_based_proposal(model, theta_f, ptcls[p], obs_i)
                        ## HACK: ([2] incorporates prior)
                        if exp(xf.log_like[2] - ptcls[p].log_like[2]) > rand()
                            mtd_gx[p] = exp(xf.log_like[3]) # - RETRIEVE MOST RECENT MARGINAL
                            ptcls[p] = xf
                            k_log[2] += 1
                            tj *= alpha
                        else
                            tj *= 0.999
                        end
                    end
                    theta[:,p] .= ptcls[p].theta
                end
                ## RB ML update
                bme[2] += log(mlr / Statistics.mean(mtd_gx))
                w .= 1  # reset w = 1
            else
                ## standard ML update
                bme[2] += log(sum(w .* gx) / sum(w))
            end # END OF sample/mutate
            # obs_min = obs_i + 1
        else
            iterate_particle!(ptcls[p], model, t, model.obs_data[obs_i])
        end
        t .= model.obs_data[obs_i].time
    end # END OF OBS LOOP
    compute_is_mu_covar!(mu, cv, theta, w)
    C_DEBUG && println(" is mcv: ", mu, cv)
    # bme[1] += log(sum(w) / outer_p)
    # bme[2] += log(Statistics.mean(w))
    bme .*= -2

    # ## MCMC sample:
    # n_iter::Int64 = (rej_sample / n_chains)
    # rejs = RejectionSample(zeros(n_chains, n_iter, length(mu)), zeros(length(mu)), zeros(length(mu), length(mu)))
    # # sample particle with weight w
    # pidx = rs_multinomial(w, n_chains)
    # k_log[1] += rej_sample
    # propd = get_prop_density(cv, propd)
    # for mc in 1:n_chains
    #     p = pidx[mc]
    #     for i in 1:n_iter
    #         ## propose new theta - independent OR conditional on current sample
    #         theta_f = ind_prop ? get_mv_param(propd, 1.0, mu) : get_mv_param(propd, tj, ptcls[p].theta)
    #         # mutate and accept with MH probability
    #         xf = model_based_proposal(model, theta_f, ptcls[p])
    #         ## HACK:
    #         if exp(xf.log_like[2] - ptcls[p].log_like[2]) > rand()
    #             ptcls[p] = xf
    #             k_log[2] += 1
    #             tj *= alpha
    #         else
    #             tj *= 0.999
    #         end
    #         # sample theta
    #         rejs.theta[mc,i,:] .= ptcls[p].theta
    #     end
    # end
    # ## compute rejection sample mu and sigma:
    # # compute_rs_mu_cv!(rejs)
    # compute_rs_mu_cv!(rejs.mu, rejs.cv, rejs.theta)
    # println(" rj mcv: ", rejs.mu, rejs.cv)
    ## return results
    output = ImportanceSample(mu, cv, theta, w, time_ns() - start_time, bme)
    println("- finished in ", Int64(round(output.run_time / C_RT_UNITS)), "s. AR = ", round(100.0 * k_log[2] / k_log[1]; sigdigits = 3), "%")
    return output
end

## MBP IBIS algorithm with dependent f/g
# - MERGE WITH MACROS ***
function run_dfg_mbp_ibis(model::HiddenMarkovModel, theta::Array{Float64, 2}, ess_rs_crit = C_DF_ESS_CRIT; n_props = 3, ind_prop = false, alpha = 1.002)
    outer_p = size(theta, 1)
    println("running MBP IBIS (DFG). n = ", outer_p)
    start_time = time_ns()
    ess_crit = ess_rs_crit * outer_p
    fn_rs = rsp_systematic

    ## initialise particles
    ptcls = Array{DFGParticle,1}(undef, outer_p)
    # theta = copy(theta_init) # GET RID? *
    for p in eachindex(ptcls)
        ic = model.fn_initial_condition()
        dfg = zeros(Int64, length(model.obs_data), length(ic))
        ptcls[p] = DFGParticle(theta[p,:], ic, copy(ic), Event[], [Distributions.logpdf(model.prior, theta[p,:]), 0.0, 0.0], dfg)
    end
    ## resampling workspace
    ptcls2 = deepcopy(ptcls)

    ## proposal distribution and scalar
    propd = Distributions.MvNormal(Matrix{Float64}(LinearAlgebra.I, size(theta,2), size(theta,2)))
    tj = 0.2

    ## initialise and run:
    w = ones(outer_p)
    mu = zeros(size(theta, 2))
    cv = zeros(size(theta, 2), size(theta, 2))
    k_log = zeros(Int64, 2)
    ## for each observation
    t = zeros(outer_p)
    model.t0_index > 0 && (t .= theta[:, model.t0_index])
    bme = zeros(2)
    gx = zeros(outer_p)
    mtd_gx = zeros(outer_p)
    for obs_i in eachindex(model.obs_data)
        if model.obs_data[obs_i].obs_id > 0
            ## for each 'outer' particle
            for p in eachindex(ptcls)
                ## compute incremental weights (i.e. run pf)
                ptcls[p].g_trans[obs_i,:] .= iterate_dfg_particle!(ptcls[p], model, t[p], model.obs_data[obs_i])
                gx[p] = exp(ptcls[p].log_like[3])
            end
            # println(" gx: ", sum(gx))
            ## COMPUTE L and update weights
            lml = log(sum(w .* gx) / sum(w))
            bme[1] += lml
            # bme[1] += log(sum(w .* gx) / sum(w))
            w .*= gx
            ##
            compute_is_mu_covar!(mu, cv, theta, w)
            # println(" y", obs_i, ": ", mu)
            ## resample and mutate if criteria satisfied:
            essv = compute_ess(w)
            # t = model.obs_data[obs_i].time
            if (essv < ess_crit)
                ## resample and swap
                propd = get_prop_density(cv, propd)
                nidx = rs_systematic(w)
                # println(obs_i, " - ", length(nidx), " - ", length(mtd_gx))
                mtd_gx .= gx[nidx]
                for p in eachindex(ptcls)
                    ptcls2[p] = deepcopy(ptcls[nidx[p]])
                end
                mlr = Statistics.mean(gx[nidx]) * exp(lml)
                ptcls, ptcls2 = ptcls2, ptcls
                # mutate:
                k_log[1] += outer_p * n_props
                for p in eachindex(ptcls)
                    for mki in 1:n_props
                        ## propose new theta - independent OR conditional on current sample (rec'd)
                        theta_f = ind_prop ? get_mv_param(propd, 1.0, mu) : get_mv_param(propd, tj, ptcls[p].theta)
                        xf = partial_dfg_model_based_proposal(model, theta_f, ptcls[p], obs_i)
                        ## HACK:
                        # if exp(sum(xf.log_like) - sum(ptcls[p].log_like)) > rand()
                        if exp(xf.log_like[2] - ptcls[p].log_like[2]) > rand()
                            mtd_gx[p] = exp(xf.log_like[3])
                            ptcls[p] = xf
                            k_log[2] += 1
                            tj *= alpha
                        else
                            tj *= 0.999
                        end
                    end
                    theta[p,:] .= ptcls[p].theta
                end
                ## RB ML update
                bme[2] += log(mlr / Statistics.mean(mtd_gx))
                w .= 1  # reset w = 1
            else
                ## standard ML update
                bme[2] += log(sum(w .* gx) / sum(w))
            end # END OF sample/mutate
            # obs_min = obs_i + 1
        else
            ptcls[p].g_trans[obs_i,:] .= iterate_dfg_particle!(ptcls[p], model, t, model.obs_data[obs_i])
        end
        t .= model.obs_data[obs_i].time
    end # END OF OBS LOOP
    compute_is_mu_covar!(mu, cv, theta, w)
    println(" is mcv: ", mu, cv)
    bme .*= -2
    ## return results
    println("- finished in ", Int64(round(output.run_time / C_RT_UNITS)), "s. AR = ", round(100.0 * k_log[2] / k_log[1]; sigdigits = 3), "%")
    return ImportanceSample(mu, cv, theta, w, time_ns() - start_time, bme)#, rejs
end
