# ModiaMath

ModiaMath provides a simulation engine and other mathematical utilities for packages 
[Modia](Modia-url) and [Modia3D](Modia3D-url)
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

see [ModiaMath installation](ModiaMathInstallation-url).


## Documentation

- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*


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
  plot(result, [(:phi, :w) :a])
```

[PendulumPlot](PendulumPlot-url)


### To run examples:
```julia
  include("$(ModiaMath.path)/examples/Simulate_Pendulum.jl")
  include("$(ModiaMath.path)/examples/Simulate_FreeBodyRotation.jl")
  include("$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_SimpleStateEvents.jl")
  include("$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_BouncingBall.jl")
```

### To run tests:
```julia
  include("$(ModiaMath.path)/test/runtests.jl")
```

## Status

The package has been tested with Julia `0.6` on Windows 7.
Version number is 0.2.0-beta.1 and functionality and robustness is planned to be improved for the 1.0 version,
see [Plans for Version 1.0](ModiaMathPlans-url).


## Issues and Contributions

Contributions are welcome, as are feature requests and suggestions.
Please open an [issue][issues-url] in this case and also if you encounter problems.


## Main Developer
[Martin Otter](https://rmc.dlr.de/sr/de/staff/martin.otter/), 
[DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)

License: MIT (expat)



[Modia-url]:   https://github.com/ModiaSim/Modia.jl
[Modia3D-url]: https://github.com/ModiaSim/Modia3D.jl
[ModiaMathInstallation-url]: https://ModiaSim.github.io/ModiaMath.jl/latest/index.html#Installation-1
[PendulumPlot-url]: https://ModiaSim.github.io/ModiaMath.jl/resources/images/pendulumPlot.png
[ModiaMathPlans-url]: https://ModiaSim.github.io/ModiaMath.jl/latest/man/Plans.html

[docs-latest-url]: https://ModiaSim.github.io/ModiaMath.jl/latest/
[docs-stable-url]: https://ModiaSim.github.io/ModiaMath.jl/stable/

[issues-url]: https://github.com/ModiaSim/ModiaMath.jl/issues
