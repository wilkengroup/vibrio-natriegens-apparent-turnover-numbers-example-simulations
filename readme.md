# Example scripts
This repo contains scripts, models, and data that can be used to reproduce the main findings of the paper:

**Multi-omic data fusion reveals the in vivo apparent turnover numbers of Vibrio natriegens at the genome scale**

## Load environment
This work exclusively used the [Julia Language](https://julialang.org/) and its
package ecosystem. Thus, to use this repo you must first install the language. Thereafter, you
must activate the environment of this repo to run the simulations (the file `Project.toml` lists the relevant information). This step uses Julia's
package manager to instantiate the correction versions of all the packages used
in the paper.

```
using Pkg
Pkg.activate(".") # activate the current repo directory is the working directory
Pkg.instantiate() # loads all packages using the Project.toml file
```

Note: this project used a commercial optimisation tool, viz.
[Gurobi](https://www.gurobi.com/). Licenses for academic use are free, but any
QP optimization tool can be used instead of `Gurobi` in case you want to try
another approach. See [here](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) for more options.

## Context specific models

## Parameter optimization example

## Simulations using optimal parameter sets


