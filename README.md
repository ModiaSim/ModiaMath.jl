# ModiaMath

[![Travis](https://travis-ci.org/ModiaSim/ModiaMath.jl.svg?branch=master)](https://travis-ci.org/ModiaSim/ModiaMath.jl)
[![AppVoyer](https://ci.appveyor.com/api/projects/status/github/ModiaSim/ModiaMath.jl?svg=true)](https://ci.appveyor.com/project/MartinOtter/modiamath-jl)
[![codecov.io](http://codecov.io/github/ModiaSim/ModiaMath.jl/coverage.svg?branch=master)](http://codecov.io/github/ModiaSim/ModiaMath.jl?branch=master)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://modiasim.github.io/ModiaMath.jl/stable/)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://modiasim.github.io/ModiaMath.jl/latest/)
[![The MIT License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square)](https://github.com/ModiaSim/ModiaMath.jl/blob/master/LICENSE.md)

ModiaMath provides a simulation engine and other mathematical utilities for packages
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

The package is registered in METADATA.jl and can be installed
in the following way (Julia >= 1.0 is required):

```julia
julia> ]add ModiaMath
```

ModiaMath uses [PyPlot](https://github.com/JuliaPy/PyPlot.jl) for plotting.
If `PyPlot` is not available in your current Julia environment
an information message is printed and all `ModiaMath.plot(..)` calls are ignored.
In order that plot windows are displayed, you need to add `PyPlot` to your current environment
via `]add PyPlot`. Often this automatic installation fails and it is recommended to follow
instead the instructions
[Installing PyPlot in a robust way](https://github.com/ModiaSim/ModiaMath.jl/wiki/Installing-PyPlot-in-a-robust-way).


## Documentation

- [**STABLE**](https://modiasim.github.io/ModiaMath.jl/stable/) &mdash; *documentation of the last released version.*
- [**LATEST**](https://modiasim.github.io/ModiaMath.jl/latest/) &mdash; *in-development version of the documentation.*


## Use

### To define a model
(note, it is simpler and less error prone to define a model with Modia or Modia3D):

```julia
  using ModiaMath
  using StaticArrays

  @component Pendulum(;L=1.0, m=1.0, d=0.1, g=9.81) begin
     phi = RealScalar(start=pi/2, unit="rad"    , fixed=true,               numericType=ModiaMath.XD_EXP)
     w   = RealScalar(start=0.0 , unit="rad/s"  , fixed=true, integral=phi, numericType=ModiaMath.XD_EXP)
     a   = RealScalar(            unit="rad/s^2",             integral=w  , numericType=ModiaMath.DER_XD_EXP)
     r   = RealSVector{2}(        unit="m"      ,                           numericType=ModiaMath.WC)
  end;

  function ModiaMath.computeVariables!(p::Pendulum, sim::ModiaMath.SimulationState)
     L = p.L; m = p.m; d = p.d; g = p.g; phi = p.phi.value; w = p.w.value

     p.a.value = (-m*g*L*sin(phi) - d*w) / (m*L^2)

     if ModiaMath.isStoreResult(sim)
        p.r.value = @SVector [L*sin(phi), -L*cos(phi)]
     end
  end;

  simulationModel = ModiaMath.SimulationModel(Pendulum(L=0.8, m=0.5, d=0.2), stopTime=5.0);

```


### To simulate a model and plot results:

```julia
  result = ModiaMath.simulate!(simulationModel; log=true);
  ModiaMath.plot(result, [(:phi, :w) :a])
```

![PendulumPlot](https://ModiaSim.github.io/ModiaMath.jl/resources/images/pendulumPlot.svg)


### To run examples and tests
```julia
  # run examples
  import ModiaMath
  include("$(ModiaMath.path)/examples/Simulate_Pendulum.jl")         # ODE as index-0 DAE
  include("$(ModiaMath.path)/examples/Simulate_FreeBodyRotation.jl") # index-1 DAE
  include("$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_PendulumDAE.jl") # index-3 DAE
  include("$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_SimpleStateEvents.jl")
  include("$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_BouncingBall.jl")
  include("$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_IdealClutch.jl") # Dirac impulses

  # run all tests
  include("$(ModiaMath.path)/test/runtests.jl")
    
```


## Status

The ModiaMath version number is 0.5.0 and functionality and robustness is planned to be improved for the 1.0 version,
see [Plans for ModiaMath version 1.0](https://ModiaSim.github.io/ModiaMath.jl/latest/man/Plans.html).


## Issues and Contributions

Contributions are welcome, as are feature requests and suggestions.
Please open an [issue](https://github.com/ModiaSim/ModiaMath.jl/issues) in this case and also if you encounter problems.


## Main Developer
[Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/),
[DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)

License: MIT (expat)
