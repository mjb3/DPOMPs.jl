## standard DA-MCMC proposals

## trajectory proposal
function std_mcmc_proposal(model::HiddenMarkovModel, xi::Particle) #, theta_f::Array{Float64,1}
    prop_type = rand(1:3)   # make proposal
    xf = deepcopy(xi)
    t0 = (model.t0_index == 0) ? 0.0 : xf.theta.value[model.t0_index]
    if prop_type == 1

    else
        ## insert / delete
        if prop_type == 1
        end
    end

    # return (Particle, prop_lk)
end


## new standard proposal function
# NEED TO BENCHMARK AGAINST OLD ***
function standard_proposal(model::PrivateDiscuitModel, xi::MarkovState, xf_parameters::ParameterProposal)
    ## choose proposal type
    prop_type = rand(1:3)
    # trajectory proposal
    xf_trajectory = deepcopy(xi.trajectory)
    t0 = (model.t0_index == 0) ? 0.0 : xf_parameters.value[model.t0_index]
    if prop_type == 3
        ## move
        length(xi.trajectory.time) == 0 && (return MarkovState(xf_parameters, xi.trajectory, NULL_LOG_LIKE, DF_PROP_LIKE, prop_type))
        # - IS THERE A MORE EFFICIENT WAY TO DO THIS? I.E. ROTATE using circshift or something?
        # choose event and define new one
        evt_i = rand(1:length(xi.trajectory.time))
        # evt_tm = (rand() * (model.obs_data.time[end] - t0)) + t0 #, xi.trajectory.event_type[evt_i])
        evt_tp = xi.trajectory.event_type[evt_i]
        # remove old one
        splice!(xf_trajectory.time, evt_i)
        splice!(xf_trajectory.event_type, evt_i)
        # add new one at random time
        add_event!(xf_trajectory, evt_tp, (rand() * (model.obs_data.time[end] - t0)) + t0)
        # compute ln g(x)
        prop_lk = 1.0
    else
        ## insert / delete
        if prop_type == 1
            ## choose type:
            tp = rand(1:size(model.m_transition, 1))
            ## insert at randomly chosen time
            add_event!(xf_trajectory, tp, (rand() * (model.obs_data.time[end] - t0)) + t0)
            ## compute ln g(x)
            prop_lk = (size(model.m_transition, 1) * (model.obs_data.time[end] - t0)) / length(xf_trajectory.time)
        else
            ## delete
            # println(" deleting... tp:", tp, " - ec: ", ec)
            length(xi.trajectory.time) == 0 && (return MarkovState(xi.parameters, xf_trajectory, NULL_LOG_LIKE, DF_PROP_LIKE, prop_type))
            # choose event index (repeat if != tp)
            evt_i = rand(1:length(xi.trajectory.time))
            # remove
            splice!(xf_trajectory.time, evt_i)
            splice!(xf_trajectory.event_type, evt_i)
            # compute ln g(x)
            prop_lk = length(xi.trajectory.time) / ((model.obs_data.time[end] - t0) * size(model.m_transition, 1))
        end # end of insert/delete
    end
    ## evaluate full likelihood for trajectory proposal and return
    return MarkovState(xi.parameters, xf_trajectory, compute_full_log_like(model, xi.parameters.value, xf_trajectory), prop_lk, prop_type)
end # end of std proposal function
