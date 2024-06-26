# son-clustering-experiments
Code to generate the figures in Dunlap and Mourrat, "Sum-of-norms clustering does not separate nearby balls"

To generate the figures, run `julia scripts/generate-plots.jl`

# More general information

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> son-clustering-experiments

It is authored by Alexander Dunlap and Jean-Christophe Mourrat.

To (locally) reproduce this project, do the following:

0. Download this code base.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "son-clustering-experiments"
```
which auto-activate the project and enable local path handling from DrWatson.
