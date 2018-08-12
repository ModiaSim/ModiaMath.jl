# ModiaMath.jl Documentation

[ModiaMath](https://github.com/ModiaSim/ModiaMath.jl) provides a simulation engine and other mathematical utilities for packages 
[Modia](https://github.com/ModiaSim/Modia.jl) and [Modia3D](https://github.com/ModiaSim/Modia3D.jl)
that are used to model physical systems such as electrical circuits, robots, or vehicles.
The recommended way is to use ModiaMath via Modia or Modia3D.
However, ModiaMath is self-contained and can be also used without Modia/Modia3D.

The central part of ModiaMath is a simulation engine to solve 
*implicit index one Differential Algebraic Equations (DAEs)*
with and without *time and state events*. The theory is partially described in 
[(Otter/Elmqvist, 2017)](http://www.ep.liu.se/ecp/132/064/ecp17132565.pdf).
In particular it is shown, that a large class of DAEs can be transformed *automatically* to this
form (including multibody systems with kinematic loops). As integrator currently 
IDA of the [Sundials integrator suite](https://computation.llnl.gov/projects/sundials)
is used (via the Julia [Sundials interface package](https://github.com/JuliaDiffEq/Sundials.jl)). 
It is planned to adapt ModiaMath to Julia package
[DifferentialEquations](https://github.com/JuliaDiffEq/DifferentialEquations.jl)
and use IDA and other appropriate integrators via this package in the future.

Additionally, ModiaMath provides functions to perform plotting in a convenient way,
to generate and use rotation matrices and quaternions for kinematic transformations in 3D, 
and to provide an infrastructure for DAE variables as needed by Modia3D.


## Installation

ModiMath is registered in METADATA.jl and can be installed with Pkg.add.

```julia
julia> Pkg.add("ModiaMath")
```

A higher level [`ModiaMath.plot`](@ref) function is provided in ModiaMath to visualize the time series
of simulation results in a convenient way. Other Julia plotting packages can also be used, but all the
tests and examples in ModiaMath, Modia and Modia3D use [`ModiaMath.plot`](@ref)). 
It is based on the Julia interface package [`PyPlot`](https://github.com/JuliaPy/PyPlot.jl) which
uses the [Matplotlib](http://matplotlib.org/) plotting library from Python. `PyPlot` need to be
installed, if the examples and tests of ModiaMath shall be executed. Installing `PyPlot` by just 
using the Julia package manager often fails. The following installation order is recommended:

1. Install a Python 3.x distribution that contains Matplotlib.\
   Recommended: [Anaconda distribution](https://www.anaconda.com/download/).\
   Advantage: very robust; disadvantage: > 3 GByte memory needed;\
   ModiaMath is based on the Python 3.x version of Matplotlib where some keywords
   are different to the Python 2.x version.
2. Include the path to the python executable in your `HOME/.juliarc.jl` file:\
    `ENV["PYTHON"] = joinpath("....", "Anaconda3", "python.exe")`
3. Start Julia, give the command `ENV["PYTHON"]` in the REPL, and check whether the path
   is correct (if you made a typo in the `.juliarc.jl` file, Julia might use another
   Python executable and PyPlot might crash Julia).
4. If you have used a different Python installation before, execute the command
   `Pkg.build["PyCall"]`, exit Julia and start Julia again.
5. Install PyPlot via `Pkg.add("PyPlot")`

