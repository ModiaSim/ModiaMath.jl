# Overview

## Mathematical Description

ModiaMath provides a simulation engine and other mathematical utilities 
to solve initial value problems of the following **implicit index 1 DAE
with events** (containing an ODE or a semi-explicit index 1 DAE as special cases;
x(t), z(t) are real vectors):

```math
\begin{align}
 z &= f_z(x, t) \\
 0 &= f_d(\dot{x}, x, t, z_i > 0) \\
 0 &= f_c(x, t, z_i > 0) \\
 J &= \left[ \frac{\partial f_d}{\partial \dot{x}};  
             \frac{\partial f_c}{\partial x} \right] \; \text{is regular}
\end{align}
```

where

```math
\lim_{\epsilon \rightarrow 0} x(t_0 - \epsilon) = x_0^{-}
```

is given. Equations ``z=z(t)`` are zero-crossing functions. Whenever a ``z_i(t)`` crosses zero 
the integration is halted, functions ``f_d, f_c`` might be changed and 
afterwards integration is restarted. If the Jacobian ``J`` is **regular**,
the DAE has an index 1 (= by differentiating ``f_c`` once, the system can be transformed
to an ODE).

Given a DAE model, ModiaMath assumes that ``J`` is **regular** for all time instants.
If this condition is violated, initialization and simulation will usually fail and an error message of the
form *"Solver does not converge"* might appear. Note, ModiaMath does not check this condition and can therefore
not provide better diagnostics in such cases.

Initial conditions ``x_0^{-}`` must be provided before simulation can start. They need 
**not** to fulfil the constraint equations, so ``f_c (x_0^{-},t_0 )≠0`` is allowed.
If this is the case, initialization will simulate for an infinitesimal small time instant 
so that ``x_0^{-}`` changes discontinuously to ``x_0^{+}`` with ``f_c (x_0^{+},t_0 )=0``. 
Note, ``\dot{x}`` is a Dirac impulse if ``x`` changes discontinuously at initialization.

As shown in xxx, every DAE can be transformed to the form above, at least in principal.
In yyy algorithms are proposed to automatically transform a large class of DAEs to this
form *without solving algebraic equations and retaining the sparsity of the equations*. 
This may require to analytically differentiating equations. The algorithms of this paper are
implemented in the Julia package `Modia` which in turn uses `ModiaMath`.
In `Modia3D` the transformation to this form is built into the package itself.

Note, the above DAE could be further transformed to an ODE (``\dot{x} = f(x,t)``), but
then the evaluation of function ``f(x,t)`` might require to solve local linear and/or
nonlinear equation systems. Furthermore, there are systems 
(for example ModiaMath/examples/Simulate_FreeBodyRotation.jl) where the
ODE states ``x`` need to be dynamically changed during simulation.

It is highly recommended to use `Modia` or `Modia3D` for simulating DAEs because this is
much simpler and less error prone as when utilizing ModiaMath directly.
However, ModiaMath can be also used without Modia or Modia3D. In this case, basically
one Julia function with the following interface has to be provided
(and in this function specific utility functions can be called)

```julia
getModelResidues(m::AbstractSimulationModel, t::Float64, x::Vector{Float64},  
                 derx::Vector{Float64}, r::Vector{Float64}
```

where `r` is the vector of residues (``r = \left[ f_d; f_c \right]``). Given the 
simulation model `m` (= a mutable struct), the actual time instant `t`, the DAE variables
`x(t)` and their derivatives `derx(t)`, the function has to compute the residue vector `r(t)`.
In directory ModiaMath/examples/withoutMacros_withoutVariables several examples are present
that are based on this interface.

In order to simplify the definition of direct ModiaMath models (to evaluate and test ModiaMath 
functionality), the macro [`@component`](@ref) has been introduced.
The examples in directory ModiaMath/examples/xxx.jl use this model definition.
In directory ModiaMath/examples/withoutMacros/xxx.jl the same examples are present,
however, the macro has been manually expanded (to show and test the result of the macro).
The [`@component`](@ref) does not yet support events. If events are present in a model,
the model has to be defined as shown in the examples of directory
ModiaMath/examples/withoutMacros_withoutVariables.


## Getting Started

You can just past the following code into the Julia REPL.


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

This results in:

```@raw html
<img src="../../resources/images/pendulumPlot.svg">
```



### To run examples:

```julia
  include("$(ModiaMath.path)/examples/Simulate_Pendulum.jl")
  include("$(ModiaMath.path)/examples/Simulate_FreeBodyRotation.jl")
  include("$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_SimpleStateEvents.jl")
  include("$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_BouncingBall.jl")
```


## Package structure

The ModiaMath package is structured in the a set of sub-modules. The most important are:

- [`ModiaMath.SimulationEngine`](@ref)\
  The engine to simulate implicit index 1 DAEs with events.

- [`ModiaMath.DAE`](@ref)\
  Interface between the SimulationEngine and the index 1 DAE model
  (e.g. initialization and event iteration is performed here).

- [`ModiaMath.Result`](@ref)\
  The `plot` function of this module allows to plot the result data of the simulation engine
  by giving the signal names. With tuples and/or vectors/matrices of signal names, the window
  layout of the figures is defined. The legends/labels of the plots are automatically constructed by
  the signal names and their units.

- [`ModiaMath.Variables`](@ref)\
  Provides Variable types to define properties of the variables on a higher level and copy
  automatically the interface vectors from the integrator into the variables and vice versa.

- [`ModiaMath.Frames`](@ref)\
  Functions that generate and operate on frames, that is coordinate systems in 3D.
  The orientation of a frame is described either with a 3x3 rotation matrix or with a 
  quaternion vector. This module is currently mainly used from Modia3D, but the functionality
  is useful for all 3D programs.




