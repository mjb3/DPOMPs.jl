# Package manual
```@contents
Pages = ["manual.md"]
Depth = 3
```

## Index
```@index
```

## Types

### Model types
```@docs
DPOMPModel
Particle
Event
Observation
```

### Results
```@docs
SimResults
ImportanceSample
RejectionSample
MCMCSample
ARQMCMCSample
```

## Functions

### Models
```@docs
generate_model
generate_custom_model
```

### Simulation
```@docs
gillespie_sim
```

### Inference

```@docs
run_mcmc_analysis
run_mbp_ibis_analysis
run_smc2_analysis
run_arq_mcmc_analysis
```

### Utilities
```@docs
get_observations
tabulate_results
print_results
```

### Visualisation

```@docs
plot_trajectory
plot_parameter_trace
plot_parameter_marginal
plot_parameter_heatmap
```

### Advanced features.
This section covers functionality for customising the algorithms themselves.

```@docs
run_custom_mcmc_analysis
generate_custom_particle
```
