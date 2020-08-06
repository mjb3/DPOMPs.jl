#### common functions ####

## choose event type
function choose_event(cum_rates::Array{Float64,1})
    etc = rand() * cum_rates[end]
    for i in 1:(length(cum_rates) - 1)
        cum_rates[i] > etc && return i
    end
    return length(cum_rates)
end

## Gaussian mv parameter proposal
function get_mv_param(propd::Distributions.MvNormal, sclr::Float64, theta_i::Array{Float64, 1})
    output = rand(propd)
    output .*= sclr
    output .+= theta_i
    return output
end

## compute mean and covar matrix for a rejection sample
# function compute_rs_mu_cv!(mu::Array{Float64,1}, cv::Array{Float64,2}, theta::Array{Float64,3})
#     # x.mu .= zeros(size(samples,1))
#     for p in eachindex(mu)
#         mu[p] = Statistics.mean(theta[:,:,p])
#     end
#     d::Int64 = length(theta) / length(mu)
#     cv .= Statistics.cov(reshape(theta, d, length(mu)))
# end

## compute mean and covar matrix for a rejection sample
function handle_rej_samples(theta::Array{Float64,3}, ap::Int64 = 0)
    # x.mu .= zeros(size(samples,1))
    output = RejectionSample(theta, zeros(size(theta,1)), zeros(size(theta,1),size(theta,1)))
    for p in 1:size(theta,1)
         output.mu[p] = Statistics.mean(theta[p,(ap+1):size(theta,2),:])
    end
    d::Int64 = size(theta,3) * (size(theta,2) - ap)
    output.cv .= Statistics.cov(transpose(reshape(theta[:,(ap+1):size(theta,2),:], size(theta,1), d)))
    return output
end
