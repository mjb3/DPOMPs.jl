# Introduction

**DPOMPs.jl** - *Bayesian inference for Discrete-state-space Partially Observed Markov Processes in Julia.*

!!! note

    Please note that this package is still in development.

## What are DPOMP models?
**Discrete-state-space (DSS)** models are used throughout ecology and other scientific domains to represent systems comprised of interacting components (e.g. people or molecules.)

A well-known example is the **Kermack-McKendrick susceptible-infectious-susceptible (SIR)** model:
```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/sir.png" alt="SIR model" style="height: 80px;"/>
```

In applied situations (e.g. like a scientific study) such systems are often difficult to directly observe, and so they are referred to in context as being **Partially Observed**.

The dynamics (how the system changes over time) of the **SIR**, and other **DSS** models, can be represented in continuous time by [a set of coupled] **Markov Processes**. Specifically, we can define a probability density, or likelihood function, that governs the time-evolution of the system under study.

Furthermore, given some [partially complete] scientific data, the combined concepts yield a paradigm for [in this case, Bayesian] statistical inference based on a general class of model:  **Discrete-state-space Partially Observed Markov Processes, or DPOMPs.**

In other words, DPOMP models and methods allow us to learn about the system of interest (e.g. an ecosystem, pandemic, chemical reaction, and so on,) even in cases where the available data is only partially complete.

## Package features

**DPOMPs.jl** is a package for:

* Bayesian parameter inference, and
* Simulation of,
* *Discrete*-state-space *Partially Observed Markov Processes*, in Julia.
* It also includes automated tools for things like convergence diagnosis, model assessment and analysis.

## Installation
See the package [source code repository](https://github.com/mjb3/DPOMPs.jl) for instructions on how to install the package.

## Overview

```@contents
Pages = [
    "models.md",
    "examples.md",
    "manual.md",
]
Depth = 2
```
