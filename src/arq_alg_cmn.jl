### common algorithm stuffs

## 'prior'
function validate_theta(theta::Array{Int64, 1}, max_idx::Int64)
    for i in eachindex(theta)
        (theta[i] < 1 || theta[i] > max_idx) && (return false)
    end
    return true
end

## get theta coords from index
function get_theta_val(model::LikelihoodModel, theta::Array{Int64, 1}, jitter = model.jitter)
    output = zeros(Float64, length(theta))
    for i in eachindex(theta)
        gap = (model.grid_range[i, 2] - model.grid_range[i, 1]) / model.sample_resolution
        output[i] = model.grid_range[i, 1] + ((theta[i] - 0.5) * gap)
        jitter > 0.0 && (output[i] += (((rand() * 2) - 1) * jitter * gap))
    end
    return output
end

## propose new theta coords
function get_theta_f(theta_i::Array{Int64, 1}, j_w::StatsBase.ProbabilityWeights, max_dist::Int64, min_dist::Int64)
    output = zeros(Int64, length(theta_i))
    d = (min_dist == max_dist) ? max_dist : rand(min_dist:max_dist)
    while sum(abs.(output)) != d# determine
        p::Int64 = StatsBase.sample(j_w)
        output[p] += Random.bitrand()[1] ? 1 : -1
        # sum(abs.(output)) == d && break
    end
    output .+= theta_i          # move
    return output
end

## inner alg constants:
const Q_JUMP = 0.1              # initial = Q_JUMP * grid res * NP
const Q_J_MIN = 2
const N_ADAPT_PERIODS = 100     # adaptive mh mcmc (parameterise?)

## adapts jump weights
function adapt_jw!(j_w::StatsBase.ProbabilityWeights, lar_j::Int64, j::Int64, mc_accepted::BitArray{1}, a_h::Int64, i::Int64, tgt_ar::Float64, mc_idx::Array{Int64,2})
    if (j == Q_J_MIN && sum(mc_accepted[(i + 1 - a_h):i]) == 0)
        C_DEBUG && print(" *LAR*")
        j = lar_j
    else
        # adjust max jump size based on acceptance rate
        j += ((sum(mc_accepted[(i + 1 - a_h):i]) / a_h) > tgt_ar ? 1 : -1)
        j = max(j, Q_J_MIN)
    end
    ## tune var
    sd = Statistics.std(mc_idx[:, 1:i], dims = 2)[:,1]
    if sum(sd) == 0.0       # adjust for zeros
        C_DEBUG && print(" *ZAR*")
        sd .= 1.0
    else
        msd = minimum(sd[sd .> 0.0])
        for s in eachindex(sd)
            sd[s] == 0.0 && (sd[s] = msd)
        end
    end
    j_w .= sd               # update weights and return j
    return j
end

## compute mean and covar matrix for a single chain
# pass mean?
function compute_chain_mean_covar(samples::Array{Float64, 3}, mc::Int64, adapt_period::Int64, steps::Int64)
    C_DEBUG && println(" - SS := ", size(samples))
    adapted = (adapt_period + 1):steps
    mc_bar = zeros(size(samples, 1))
    for i in eachindex(mc_bar)
        mc_bar[i] = Statistics.mean(samples[i, adapted, mc])
    end
    scv = Statistics.cov(transpose(samples[:, adapted, mc]))
    C_DEBUG && print(" - scv := ", scv)
    return (mc_bar, scv)
end

## initialise inner ARQ MCMC
macro init_inner_mcmc()
      esc(quote
      ## adaption interval
      a_h::Int64 = max(steps / N_ADAPT_PERIODS, 100)
      ## adaptive stuff
      j::Int64 = round(Q_JUMP * model.sample_resolution * length(theta_i))
      j_w = StatsBase.ProbabilityWeights(ones(length(theta_i)))
      ## declare results
      mc_idx = Array{Int64, 2}(undef, length(theta_i), steps)
      # mc = Array{Float64, 2}(undef, steps, length(theta_init))
      mc_log_like = Array{Float64,1}(undef, steps)
      mc_accepted = falses(steps)
      mc_time = zeros(UInt64, steps)
      # DEBUG (TO BE REMOVED):
      mcf = Array{Float64, 2}(undef, length(theta_i), steps)
      ## theta index (i.e. grid key)
      ## first sample
      mc_idx[:,1] .= theta_i
      mc_accepted[1] = true
      mc_fx = zeros(Int16, steps)   # process run

      ## estimate x0
      xi = get_grid_point!(grid, theta_i, model, true)
      ## write first sample and run the Markov chain:
      samples[:,1,mc] .= xi.result.sample
      mcf[:,1] .= xi.result.sample
      mc_log_like[1] = xi.result.log_likelihood
      mc_fx[1] = xi.process_run
      mc_time[1] = time_ns()
      end)
end
