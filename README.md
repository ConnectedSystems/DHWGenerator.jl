# DHWGenerator.jl

Generator of synthetic DHW trajectories.

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
