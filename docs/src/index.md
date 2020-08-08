# Introduction

**DPOMPs.jl** - *Bayesian parameter inference for Discrete-state-space Partially Observed Markov Processes in Julia.*

!!! note

    Please note that this package is still in development.

## What are DPOMP models?
**Discrete-state-space (DSS)** models are used throughout ecology and other domains to represent systems of interacting components (e.g. people or molecules.)

A well-known example is the Kermack-McKendrick SIR model:
```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/sir.png" alt="SIR model" style="height: 60px;"/>
```

In applied scientific situations, such systems are often difficult to directly observe, and so they are referred to in context as **Partially Observed**.

Finally, the dynamics (how the system state changes) of the **SIR** and other **DSS** models can be represented in continuous time by [a set of coupled] **Markov Processes**. Specifically, we can define a probability density or likelihood function.

This is useful outwith the availability of scientific data, because we can construct statistical sampling schemes based on that identity (i.e. the probability density function.) In plain English - it allows us to 'simulate' the model, and thus gain an intuitive understand of those dynamics.

Given some (partially complete) data however, these concepts yield a paradigm for (in this case, Bayesian) statistical inference based on a general class of model:  **Discrete-state-space Partially Observed Markov Processes, or DPOMPs.**

## Package features

**DPOMPs.jl** is a package for:

* Bayesian parameter inference, and
* Simulation of,
* Discrete-state-space Partially Observed Markov Processes, in Julia.
* It also includes automated tools for convergence diagnosis and analysis.

## Installation
See the package [source code repository][https://github.com/mjb3/DPOMPs.jl] for instructions.

## Overview

```@contents
Pages = [
    "models.md",
    "examples.md",
    "manual.md",
]
Depth = 2
```
