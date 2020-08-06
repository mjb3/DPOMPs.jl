### standard ARQ MCMC algorithm

## standard grid request
function get_grid_point!(grid, theta_i::Array{Int64, 1}, model::LikelihoodModel, burn_in::Bool)  #, herd_size_dist::EmpiricalDist, log_irw_y::Array{Float64, 2}
    Q_BI_SAMPLE_LIM = 1
    ## check hashmap
    exists = haskey(grid, theta_i)
    if exists
        # return log value if visited has hit max
        x = grid[theta_i]
        visited = x.visited
        sampled = x.sampled
        # visited > (burn_in ? Q_BI_SAMPLE_LIM : model.sample_limit) && (return GridRequest(x, false))
        theta_val = x.sample
    else
        visited = 0
        sampled = 0
        theta_val = get_theta_val(model, theta_i)
    end
    ## VISITED: COMPUTED
    if visited < (burn_in ? Q_BI_SAMPLE_LIM : model.sample_limit)
        ## update density
        log_like = model.pdf(theta_val)
        exists && (log_like = log(exp(x.log_likelihood) + ((exp(log_like) - exp(x.log_likelihood)) / visited)))
        visited += 1
        comp = true
    else
        log_like = x.log_likelihood
        comp = false
    end
    burn_in || (sampled += 1)
    ## update hashmap
    output = GridPoint(theta_val, log_like, visited, sampled)
    grid[theta_i] = output
    ## return updated sample
    return GridRequest(output, comp)
end

## run standard inner MCMC procedure and return results
function arq_met_hastings!(samples::Array{Float64,3}, mc::Int64, grid::Dict, model::LikelihoodModel, steps::Int64, adapt_period::Int64, theta_i::Array{Int64, 1}, tgt_ar::Float64) #, prop_type::Int64
    C_LAR_J_MP = 0.2
    lar_j::Int64 = round(C_LAR_J_MP * model.grid_resolution * length(theta_i))
    ## initialise inner MCMC
    @init_inner_mcmc
    C_DEBUG && print("- mc", mc, " initialised ")
    for i in 2:steps
        ## propose new theta
        theta_f = get_theta_f(theta_i, j_w, j, 1)
        ## validate
        if validate_theta(theta_f, model.grid_resolution)
            ## get log likelihood
            xf = get_grid_point!(grid, theta_f, model, i < a_h) # limit sample (n=1) for first interval only
            mcf[:,i] .= xf.result.sample
            mc_fx[i] = xf.process_run ? 1 : 0
            ###
            mc_log_like[i] = xf.result.log_likelihood
            ## mh step
            mh_prob = exp(mc_log_like[i] - mc_log_like[i - 1])
            # accept or reject
            if mh_prob > 1.0
                mc_accepted[i] = true
            else
                mh_prob > rand() && (mc_accepted[i] = true)
            end
        else
            ## reject automatically
            mcf[:,i] .= samples[:,i-1,mc] # temp fix: backfill with theta i-1 TBR
        end
        ## acceptance handling
        if mc_accepted[i]
            samples[:,i,mc] .= xf.result.sample
            mc_idx[:,i] .= theta_f
            theta_i = theta_f
        else
            ## rejected
            samples[:,i,mc] .= samples[:,i-1,mc]
            mc_idx[:,i] .= mc_idx[:,i - 1]
            ## TRIGGER REFRESH (XI) if rejected N times
            Q_REJECT_TRIGGER = 100
            if sum(mc_accepted[max(1, i - Q_REJECT_TRIGGER):i]) == 0
                xi = get_grid_point!(grid, theta_i, model, false)
                xi.process_run && (mc_fx[i] += 1)
                mc_log_like[i] = xi.result.log_likelihood
            else
                mc_log_like[i] = mc_log_like[i - 1]
            end
        end ## end of acceptance handling
        mc_time[i] = time_ns() - mc_time[1]
        ## ADAPTATION
        i % a_h == 0 && (j = adapt_jw!(j_w, lar_j, j, mc_accepted, a_h, i, tgt_ar, mc_idx))
        ## end of adaption period
    end ## end of Markov chain for loop
    if C_DEBUG
        print("- mc", mc, " processing -> ")
        ## compute mean/var and return results
        mc_m_cv = compute_chain_mean_covar(samples, mc, adapt_period, steps)
        C_DEBUG && print("- mcv := ", mc_m_cv, " ")
    end
    ## compute SMC runs and return results
    # return MCMCResults(sum(mc_fx), mc_idx, mc_accepted, mc_fx, mc_m_cv[1], mc_m_cv[2], mcf, mc_log_like, mc_time)
    return (sum(mc_fx), sum(mc_accepted) / steps, sum(mc_accepted[(adapt_period+1):steps]) / (steps - adapt_period))
end
