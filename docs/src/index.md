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
The latest released version (0.3.1) is the last one with support for Julia >= 0.6.
Trunk and later versions support Julia >=1.0.

```julia
# Julia >= 0.6:
julia> Pkg.add("ModiaMath")

# alternatively in Julia >= 0.7:
julia> ]add ModiaMath
```

ModiaMath uses PyPlot for plotting.
If `PyPlot` is not available in your current Julia environment
an information message is printed and all `ModiaMath.plot(..)` calls are ignored.
In order that plot windows are displayed, you need to add `PyPlot` to your current environment
via `Pkg.add("PyPlot")`. Often this automatic installation fails and it is recommended to follow
instead the instructions
[Installing PyPlot in a robust way](https://github.com/ModiaSim/ModiaMath.jl/wiki/Installing-PyPlot-in-a-robust-way).


## Release Notes


### Version 0.4.0

- All Julia 0.6 code removed.
- ModiaMath.plot supports variables that can be computed from a vector of struct result
  (for details, see docu of plot)


### Version 0.3.1

- Initialization issue corrected
- Use PyPlot.pygui(true) when importing ModiaMath in order that separate plots windows are used and no inline plots, independent of the used environment (especially made so that plots in vs-code are not inlined in the vs-code GUI).
- ModiaMath.plot supports now also Vectors of structs. E.g. Vector{ThermodynamicState} and plot(result, "state.p") is then possible, where state is an instance of the struct and "p" is a fieldname with a scalar Number value.


### Version 0.3.0

- ModiaMath initialization and re-initialization improved:
  * Initialization and re-initialization considerably changed
    (now only one nonlinear equation is solved and no longer two nonlinear equations).
  * By this change, several issues in the initialization were fixed.
  * x\_nominal also supported for ModiaMath.Variables.


### Version 0.2.6

- ModiaMath initialization and simulation improved:
  * Issue with defaultInterval fixed:
    Previously: defaultInterval = (defaultStopTime - defaultStartTime)/500.0. However, if stopTime
    was set to 10e5, then the defaultInterval was much too small. This was changed so that
    defaultInterval is either explicitly given, or it is computed from the actual values of StopTime and StartTime.
  * x\_nominal introduced in constructor of simulation state (DAE.SimulationState) and of
    Modia interface (ModiaToModiaMath.ModiaSimulationModel):
    If x\_nominal is explicitly provided in these constructor calls, it is used.
    Otherwise `x_nominal[i]` is set to `max(abs(x_start[i]), 1e-7)`.
  * x\_nominal is newly used in KINSOL to scale the unknowns and
    is used to compute the absolute tolerance for Sundials IDA.
  * A bug was corrected in the initialization function, where rScale was reported to KINSOL,
    but was also used for scaling of the residue (although this is performed in KINSOL).
  * If log=true, name-, start-, fixed-, nominal-values of the x-vector are printed before
    initialization starts.

- New function *ModiaMath.solveOneNonlinearEquation*:
  Determines the solution of one non-linear algebraic equation `y=f(u)`
  in one unknown `u` in a reliable and efficient way (using Brents algorithm).



### Version 0.2.5

- ModiaMath initialization improved:
  * KINSOL tolerance FTOL takes IDA tolerance into acount.
  * More information about the model is provided in case KINSOL fails.

- Interface to Modia improved:
  New keyword option hev (step size used for implicit Euler during initialization and at events) added to "simulate".

- ModiaMath result handling improved:
  * The result may contain single values (such as parameter J=0.1).
    When using plot(..) on such a result variable, a constant line is plotted
    between the first and the last point of the x-axis.
  * New function `ModiaMath.resultTable(result)` to transform the variable information
    in the result data structure in a DataFrames table
    (containing variable name + type + size ) that can then be printed.
  * New functions `ModiaMath.closeFigure(figure)` and
    `ModiaMath.closeAllFigures()` to close a specific figure or close all figures.
  * Changing some default options of PyPlot. In particular, (a) the labels on the x- and y-axis
    use exponential notation (e.g. 1e5), if the numbers are larger as 1e4 or smaller as 1e-3
    (PyPlot default is 1e7 and 1e-7), (b) smaller fonts and linewidth are used.
    Default options are changed via PyCall. If PyCall is not in the Julia environment, PyCall is added.
  * If a vector or matrix of subplots is defined, then the x-axis labels
    are only displayed for the last subplot row.
  * Result data is always converted to Float64 before passing it to PyPlot.
    Therefore, Julia numbers can be plotted, even if the types of the numbers are not supported by PyPlot
    (for example rational numbers can be plotted).


### Version 0.2.4

- **Non-backwards compatible change of ModiaMath.plot(..)**:
  The function was changed to only support result dictionaries of the type `Dict{AbstractString,Any}` to
  simplify implementation and maintenance and to use the identical dictionary type as in Modia.

- ModiaMath.plot(..) improved:
  * For the dictionary type used by Modia, all elements of a variable vector
    are plotted by giving the variable name, or one specific element of the variable vector is plotted by giving
    the name and the vector index.
  * A new keyword argument `maxLegend=10` introduced: If the number of entries in a legend exceeds this number,
    no legend is included in the plot. Note, the curves can still be identified by clicking in the tool bar of the
    plot window on button `Edit axis, curve ..`.
  * If the plot function is used in the wrong way (e.g. signals shall be plotted that are not in the result dictionary, or selecting a
    variable vector as x-axis), warning messages are printed and the call is ignored.

- Dependent packages updated to their newest versions.
  Especially, warnings from Unitful do no longer occur, due to the update to version 0.12.0.

- Issue with PyPlot warnings fixed:
  When ModiaMath was used in another package (e.g. in Modia or Modia3D),
  strange warning messages appeared when importing it.
  This issue was fixed with the solution sketched in
  [discourse](https://discourse.julialang.org/t/how-to-remove-pyplot-jl-from-the-distribution-of-a-package-that-uses-it/15130/10?u=martinotter).


### Version 0.2.3

- Wrong UUID of ModiaMath corrected (did not correspond to the UUID in Julias METADATA).


### Version 0.2.2

- PyPlot was removed from the REQUIRE and Project.toml files and code was added,
  so that PyPlot is automatically imported in ModiaMath if it is available in
  the current environment of the user.
  The benefit is that ModiaMath can be used, even if PyPlot is not installed.
  This is especially useful for ContinuousIntegration, because automatic
  installation of PyPlot often fails.
  Inspect the wiki page
  [Installing PyPlot in a robust way](https://github.com/ModiaSim/ModiaMath.jl/wiki/Installing-PyPlot-in-a-robust-way)
  to install PyPlot in a robust way.

- All extra packages used in examples and tests are now referenced via ModiaMath
  (for example `using ModiaMath.StaticArrays` instead of `using StaticArrays`).
  The benefit is that all examples and tests can be directly executed with `include`
  (for example: `import ModiaMath; include($(ModiaMath.path)/examples/Simulate_Pendulum.jl)`)
  provided `ModiaMath` is in the current environment. Previously, it was assumed that these
  extra packages are present in the users environment and an error occured, if this was not the case.

- New arguments `prefix` and `reuse` added to ModiaMath.plot(..) to add new plots
  (e.g. from a new simulation run) to an existing figure, without clearing the figure beforehand.


### Version 0.2.1

- Adapted to Julia 0.7 and 1.0
  (including using new package manager via Project.toml, Manifest.toml files
  and adapting the README.md and documentation files).

- Buttons for Travis CL, code coverage and docs added to README.md file.
  (this includes adaptations to .travis.yml file). Note, Appveyor CL is not
  yet activated (although appveyor.yml file is present).

- Interface call to KINSOL slightly improved to avoid a gc-crash.

- In case KINSOL fails, more run-time information is added to the error message.

- Indentation was changed consistently to 4 spaces.


### Version 0.2.0

- First public release (for Julia 0.6)
