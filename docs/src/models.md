# Models

## Guide

### Predefined models

### Customising predefined models

### Custom models from scratch

## Model directory

### Epidemiological models

#### SIR model
The canonical Kermack-McKendrick susceptible-infectious-recovered model is perhaps the best known example of state-space models within epidemiology.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/sir.png" alt="SIR model" style="height: 80px;"/>
```

#### SI model
The susceptible-infectious model is the simplest conceptual example; two states and only one type event.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/si.png" alt="SI model" style="height: 80px;"/>
```

#### SIS model
Another common derivative of the SIR model.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/sis.png" alt="SI model" style="height: 80px;"/>
```

#### SEI model
The SEI model includes an 'exposed' state, i.e. for modelling communicable diseases with *latent* non-infectious periods.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/sei.png" alt="SI model" style="height: 80px;"/>
```

#### SEIR model
Somewhat obviously, the SEIR model concept combines the SEI with the SIR.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/seir.png" alt="SI model" style="height: 80px;"/>
```

### Others

#### The Lotka-Volterra predator-prey model

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/lotka.png" alt="SIR model" style="height: 80px;"/>
```

#### Ross-MacDonald two-species Malaria model

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/img/rossmac.png" alt="SIR model" style="height: 100px;"/>
```
