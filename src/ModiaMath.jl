# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control


"""
ModiaMath - Mathematical Utilities for Modia and Modia3D

To define a model use Modia or Modia3D. You can define a model
natively in ModiaMath in the following way:

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

To simulate a model and plot results:

```julia
  result = ModiaMath.simulate!(simulationModel; log=true);
  ModiaMath.plot(result, [(:phi, :w) :a])
```

[PendulumPlot](PendulumPlot-url)


To run examples:
```julia
  include("$(ModiaMath.path)/examples/Simulate_Pendulum.jl")
  include("$(ModiaMath.path)/examples/Simulate_FreeBodyRotation.jl")
  include("$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_SimpleStateEvents.jl")
  include("$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_BouncingBall.jl")
```

To run tests:
```julia
  include("$(ModiaMath.path)/test/runtests.jl")
```

For more information, see (https://github.com/ModiaSim/ModiaMath.jl/blob/master/README.md)
"""
module ModiaMath

# Exported symbols
export @component
export RealVariable, RealScalar, RealSVector, RealSVector3
export plot


# Types
"""
    type ModiaMath.AbstractSimulationModel

Struct that is used as simulation model (has field: simulationState)
"""
abstract type AbstractSimulationModel end


"""
    type ModiaMath.AbstractComponentWithVariables

Struct that contains ModiaMath.AbstractVariables as field or as field
in a sub-struct.
"""
abstract type AbstractComponentWithVariables end

" The internal part of a component (has at least fields \"name\" and \"within\") "
abstract type AbstractComponentInternal end


"""
    type ModiaMath.AbstractVariable <: ModiMath.AbstractComponentWithVariables

A Variable used as element of the DAE model description and is 
included in the result (if no residue)
"""
abstract type AbstractVariable <: AbstractComponentWithVariables end


"""
    ModiaMath.AbstractRealVariable <: ModiaMath.AbstractVariable

A real [`ModiaMath.AbstractVariable`](@ref) (either scalar or array)
"""
abstract type AbstractRealVariable <: AbstractVariable end


"""
    const path

Absolute path of package directory of ModiaMath
"""
const path = dirname(dirname(@__FILE__))   # Absolute path of package directory
const Time = Float64   # Prepare for later Integer type of time
const Version = "0.2.5-dev from 2018-10-07 15:34"

println(" \nImporting ModiaMath version ", Version)


getVariableAndResidueValues(extraInfo::Any) = nothing    # Return a variable value and a residue value table of nonlinear solver (for error message)


# include sub-modules and make symbols available that have been exported in sub-modules
include("Utilities.jl")
using .Utilities

include("Logging.jl")
using .Logging

include(joinpath("Frames", "_module.jl"))
using .Frames

#include("SparseJacobian.jl")
#using .SparseJacobian

include(joinpath("Result", "_module.jl"))
using .Result

include(joinpath("NonlinearEquations", "_module.jl"))
const NonlinearEquationsInfo   = NonlinearEquations.KINSOL.NonlinearEquationsInfo
const solveNonlinearEquations! = NonlinearEquations.KINSOL.solveNonlinearEquations!

include(joinpath("DAE", "_module.jl"))
using .DAE

include(joinpath("SimulationEngine", "_module.jl"))
using .SimulationEngine

include(joinpath("Variables"       , "_module.jl"))
using .Variables

include("ModiaToModiaMath.jl")
using .ModiaToModiaMath



# Import packages that are used in examples and tests
# (in order that there are no requirements on the environment 
#  in which the examples and tests are executed).
import DataFrames
import StaticArrays
import Unitful

@static if VERSION >= v"0.7.0-DEV.2005"
    import LinearAlgebra
    import Test
end



end # module 
