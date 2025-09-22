# DHWGenerator.jl

Synthetic DHW generation for coral reefs.

This package provides methods to generate synthetic DHW trajectories to aid in model
testing, validation, and scenario exploration.

Chiefly, the methods can be used to create indicative time series that follow an assumed:

- Rate of warming
- Year-on-year variability (seasonal and short-term fluctuations)
- Number of extreme events

## Setup

This package is not yet registered.

If using in production, install directly from the GitHub repository:

```julia
] add https://github.com/ConnectedSystems/DHWGenerator.jl
```

## Development setup

As usual, clone project from GitHub.

Create a `sandbox` directory inside the `DHWGenerator.jl` project.

Start julia (or use the VS Code REPL)

```shell
julia --project=.
```

```julia-repl
] dev ..
```

Add plotting packages if desired:

```julia-repl
] add GLMakie
```

## Usage examples

```julia
using DHWGenerator

# Generate DHW trajectories for 20 years and 1 location
generate_dhw_trajectories(20, 1)
```

Other options can be viewed through the usual help prompt:

```julia
?generate_dhw_trajectories
```
