# DPOMPs.jl
Fast Bayesian parameter inference for Discrete-state-space Partially Observed Markov Processes in Julia

https://github.com/mjb3/DPOMPs.jl/workflows/Documentation/badge.svg
https://github.com/mjb3/DPOMPs.jl/workflows/TagBot/badge.svg

This package contains tools for Bayesian inference and simulation of DPOMP models. See the 

## Features

- Customisable finite adaptive MCMC algorithm for fast parameter inference.
- Pooley model based proposal (MBP) method for improved mixing.
- Simulation via the Gillespie direct method.
- Automated tools for convergence diagnosis and analysis.

## Installation

The package is not registered and must be added via the package manager Pkg.
From the REPL type `]` to enter the Pkg mode and run:

```
pkg> add https://github.com/mjb3/DPOMPs.jl
```

## Usage

See the [package documentation][docs] for further information and examples.

[docs]: https://mjb3.github.io/DPOMPs.jl/
