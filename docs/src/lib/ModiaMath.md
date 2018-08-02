# ModiaMath Constants and Types

The following **constants** are defined in ModiaMath:

- `const path`\
  Absolute path of ModiaMath package directory. This allows for example to run ModiaMath examples
  as `include("$(ModiaMath.path)/examples/Simulate_Pendulum.jl")`.


The following **abstract types** are defined and used in ModiaMath:

- `ModiaMath.AbstractSimulationModel`\
  Struct that is used as simulation model (has field 
  `simulationState::ModiaMath.SimulationState`). An instance of this
  type can be directly used in a [`ModiaMath.simulate!`](@ref)`(...)` call.

- `ModiaMath.AbstractComponentWithVariables`\
   Struct that contains ModiaMath.AbstractVariables as field or as field
   in a sub-struct. 

- `ModiaMath.AbstractVariable <: ModiMath.AbstractComponentWithVariables`\
   A Variable used as element of the DAE model.

- `ModiaMath.AbstractRealVariable <: AbstractVariable`\
   A real `ModiaMath.AbstractVariable` (either scalar or array).



