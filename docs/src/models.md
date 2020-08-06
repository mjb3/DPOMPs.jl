# Models

## Guide directory

### Predefined models

### Customising predefined models

### Starting from scratch

## Model directory

### Epidemiological models

#### SIR model
The canonical Kermack-McKendrick susceptible-infectious-recovered model is perhaps the best known example of state-space models within epidemiology.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/assets/mdl/sir.png" alt="SIR model" height="100"/>
```

#### SI model
The susceptible-infectious model is the simplest conceptual example; two states and only one type event.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/assets/mdl/si.png" alt="SI model" height="100"/>
```

#### SIS model
Another common derivative of the SIR model.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/assets/mdl/sis.png" alt="SI model" height="100"/>
```

#### SEI model
The SEI model includes an 'exposed' state, i.e. for modelling communicable diseases with *latent* non-infectious periods.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/assets/mdl/sei.png" alt="SI model" height="100"/>
```

#### SEIR model
Somewhat obviously, the SEIR model concept combines the SEI with the SIR.

```@raw html
<img src="https://raw.githubusercontent.com/mjb3/DPOMPs.jl/master/docs/assets/mdl/seir.png" alt="SI model" height="100"/>
```

### Others

## Custom models
