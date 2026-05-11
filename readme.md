# Example scripts
This repo contains scripts, models, and data that can be used to reproduce the main findings of the paper:

**Multi-omic data fusion reveals the in vivo apparent turnover numbers of Vibrio natriegens at the genome scale**

## Load environment
This work exclusively used the [Julia Language](https://julialang.org/) and its
package ecosystem. Thus, to use this repo you must first install the language.
Thereafter, you must activate the environment of this repo to run the
simulations. The file `Project.toml` contains the relevant information, and is
automatically used by Julia's package manager to instantiate the correct
environment. 

```
using Pkg
Pkg.activate(".") # activate the current repo directory is the working directory
Pkg.instantiate() # loads all packages using the Project.toml file
```

Note: you will need manually install the Vibrio natriegens model package like this:
```
Pkg.add("https://github.com/wilkengroup/VibrioNatriegens.git")
```
The reason for this is that the model is not a registered package.

Note: this project used a commercial optimisation tool, viz.
[Gurobi](https://www.gurobi.com/). Licenses for academic use are free, but any
QP/LP optimization solver can be used instead of `Gurobi` in case you want to
try another approach. See
[here](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) for more
options.

## Prepare data
The folder `data` contains prepared molar masses,  models, flux, and proteomic
data that can be directly integrated into the condition specific models.
However, to simplify the examples in `scripts` it is first necessary to prepare
the data for analysis. Run `scripts/prepare-data.jl` and this will generate
serialized data into  `tmp`.

## Parameter optimization example
The script `scripts/parameter-optimization.jl` contains all the code necessary
to find kapps for each condition that best explain the biological observations.  

## Simulations using optimal parameter sets
The script `scripts/run-example-ecfba.jl` runs a simulation predicting the
growth rate using `kmax` or `kapp` based kinetics.


