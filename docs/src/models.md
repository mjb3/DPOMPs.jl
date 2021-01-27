# Models

This section provides instructions for generating models (WIP).

## Guide

### Predefined models

### Customising predefined models

### Custom models from scratch

## Model directory

### Epidemiological models

#### SIR model
The canonical Kermack-McKendrick susceptible-infectious-recovered model is perhaps the best known example of state-space models used within the field of epidemiology.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/sir.png" alt="SIR model" style="height: 80px;"/>
```

```@repl
using DPOMPs
generate_model("SIR", [100, 1, 0])
```

#### SI model
The susceptible-infectious model is the simplest conceptual example of this class of model; two states and only one type of event.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/si.png" alt="SI model" style="height: 80px;"/>
```

```@repl
generate_model("SI", [100, 1])
```

#### SIS model
Another common derivative of the SIR model.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/sis.png" alt="SIS model" style="height: 90px;"/>
```

```@repl
generate_model("SIS", [100, 1])
```

#### SEI model
The SEI model includes an 'exposed' state, i.e. for modelling communicable diseases with *latent* non-infectious periods.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/sei.png" alt="SEI model" style="height: 80px;"/>
```
```@repl
generate_model("SEI", [100, 0, 1])
```

#### SEIR model
Somewhat obviously, the SEIR model concept combines the SEI with the SIR.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/seir.png" alt="SEIR model" style="height: 80px;"/>
```

```@repl
generate_model("SEIR", [100, 0, 1, 0])
```

### Others

#### The Lotka-Volterra predator-prey model

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/lotka.png" alt="Lotka model" style="height: 80px;"/>
```

```@repl
generate_model("LOTKA", [70, 70])
```

#### Ross-MacDonald two-species Malaria model

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/rossmac.png" alt="Malaria model" style="height: 160px;"/>
```

```@repl
generate_model("ROSSMAC", [100, 0, 400, 50])
```
