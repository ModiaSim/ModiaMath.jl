# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module
#   ModiaMath.DAE (ModiaMath/DAE/_module.jl)
#

# Utility functions used below
"""
    defaultVariableName(model, vcat, vindex)

Return default names for the variables (e.g. `x[1]`)
"""
defaultVariableName(model::Any, vcat::VariableCategory, vindex::Int) = vcat == Category_X ? "x[" * string(vindex) * "]" :
                                                                      (vcat == Category_W ? "w[" * string(vindex) * "]" :
                                                                                        "der(x[" * string(vindex) * "])" )

#nameVector(name, nNames::Int) = [Symbol(name, "[", i, "]") for i = 1:nNames]
#fcNameVector(fc, names::AbstractVector) = [ Symbol(fc, "(", names[i], ")") for i in eachindex(names)]

nameVector(name, nNames::Int) = [string(name, "[", i, "]") for i = 1:nNames]
fcNameVector(fc, names::AbstractVector) = [ string(fc, "(", names[i], ")") for i in eachindex(names)]

"""
    getVariableName(model, vcat, vindex, nx=0;
                    xNames   =nameVector("x", nx),
                    derxNames=fcNameVector("der", xNames),
                    wNames   =String[])

Given category `vcat` and index `vindex`of category,
return the full path name of the respective variable.
"""
function getVariableName(model::Any, vcat::VariableCategory, vindex::Int, nx=0;
                         xNames::Vector{String}=nameVector("x", nx),
                         derxNames::Vector{String}=fcNameVector("der", xNames),
                         wNames::Vector{String}=String[])
    if vcat == Category_X
        return xNames[vindex]
    elseif vcat == Category_DERX
        return derxNames[vindex]
    elseif vcat == Category_W
        return wNames[vindex]
    else
        error("Wrong argument vcat = ", vcat)
    end
end

"""
   getResultNames(model)

Return a vector of result names.
"""
function getResultNames(model::Any)
    sim::SimulationState = model.simulationState
    resultNames = ["time"; String[sim.getVariableName(model, Category_X, i)  for i = 1:sim.nx];
                          String[sim.getVariableName(sim, Category_DERX, i) for i = 1:sim.nx];
                          String[sim.getVariableName(sim, Category_W, i)    for i = 1:sim.nw] ]
    return resultNames
end


"""
    @enum StructureOfDAE
          DAE_LinearDerivativesAndConstraints
          DAE_ExplicitDerivatives
          DAE_NoSpecialStructure

Enumeration defining the structure of the DAE of the simulation model.
The following DAE structures are supported
(function `getModelResidues!(model, t, x, derx, r, w; simulation=true)` returns the residues `r`):


# ModiaMath.DAE_LinearDerivativesAndConstraints (default)

```math
\\begin{align}
     z &= f_z(x,t) \\\\
 0 = r &= \\left[ \\begin{array}{l}
                    M_d(x,t,z_i>0) \\cdot \\dot{x} + b_d(x,t,z_i>0) \\;\\; (= f_d) \\\\
                    f_c(x,t,z_i>0)
                  \\end{array} \\right] \\\\
     J &= \\left[ M_d;
                  \\frac{\\partial f_c}{\\partial x} \\right] \\; \\text{is regular (matrix is invertible)}
\\end{align}
```

where

```math
\\lim_{\\epsilon \\rightarrow 0} x(t_0 - \\epsilon) = x_0^{-}
```

Equations ``f_d`` are linear in the derivatives ``\\dot{x}``.
It is required that the Jacobian ``J`` is **regular**, that is the DAE has an index 1
(= by differentiating ``f_c`` once, the system can be transformed to an ODE).

Equations ``z=z(t)`` are zero-crossing functions. Whenever a ``z_i(t)`` crosses zero,
an event is triggered and simulation is halted. During an event, ``z_i > 0`` can
change its value. The equations above are solved with a fixed-point iteration scheme (= *event iteration*)
until ``z_i > 0`` does not change anymore. Afterwards, integration is
restarted and the Boolean variable ``z_{pos} = z_{i,ev} > 0`` keeps its value until
the next event occurs.

At an event instant some ``f_c`` equations might become ``f_d`` equations and vice versa.
The constraint equations ``f_c`` can be at any position of the residue vector `r` and at an event
instant ``f_c`` equations might become ``f_d`` equations and vice versa.
When instantiating a [`SimulationState`](@ref), the initial definition
of the constraint equations is provided with vector argument `is_constraint`.
Note, it is also possible to define time events, so triggering events at pre-defined time
instants, for example to model sampled data systems
(see [`ModiaMath.setNextEvent!`](@ref)).

Initial conditions ``x_{ev}^{-}`` must be provided before simulation can start (``x_{ev}^{-} = x_0^{-}``) or at
an event restart. They need
**not** to fulfill the constraint equations, so ``f_c (x_{ev}^{-},t_{ev} ) \\neq 0 `` is allowed.
If this is the case, initialization/re-initialization will simulate for an infinitesimal small time instant
so that ``x_{ev}^{-}`` changes discontinuously to ``x_{ev}^{+}`` with ``f_c (x_{ev}^{+},t_{ev} )=0``.
This is performed by *analytically* integrating over the initial time or the event time and
might imply to integrate over Dirac impulses (in case ``x`` is discontinuous at this time instant).
Under certain conditions a numerical approximation of the
mathematical (exact) solution is computed.

The derivative of the constraint equations ``f_c`` can be provided at
event restart. In this case ``\\frac{\\partial f_c}{\\partial x}`` is computed
with the provided ``\\dot{f}_c``, instead of computing it with finite differences
(which is numerically less reliable). In both cases ``\\dot{x}_{ev}^{+}`` is then
computed by solving a linear equation system with Jacobian ``J``.


# ModiaMath.DAE_ExplicitDerivatives

```math
\\begin{align}
     z &= f_z(x, t) \\\\
 0 = r &= - \\dot{x} + b_d(x,t,z_i>0) \\;\\; (= f_d)
\\end{align}
```
where

```math
x(t_0) = x_0
```

This is a special case of the first form, where no constraint equations are present
and the derivatives are explicit. This form is also called ODE (Ordinary Differential Equations
in state space form). It has the advantage that many ODE integration methods can be used
to solve this system.

Initialization and re-initialization is trivial, because ``x_{ev}^{+}`` is provided
as initial value or at event restart from the model and then:

```math
der(x_{ev})^{+} := f_d(0, x_{ev}^{+}, t_{ev})
```

Note, if a Dirac impulse occurs in the model, then this property has to be
handled inside the model and the result of analytically integrating over the event
instant must be returned as ``x_{ev}^{+}`` from the model.


# ModiaMath.DAE_NoSpecialStructure

```math
\\begin{align}
     z &= f_z(x,t) \\\\
 0 = r &= f(\\dot{x},x,t,z_i>0)
 \\end{align}
```

This is a general DAE. Initialization and re-initialization is performed by using an implicit Euler step.
When appropriately scaling `r`, using a step size that tends to zero and under further
assumptions, then this solution can be interpreted as analytically integrating over
the time instant. This might mean to integrate over
Dirac impulses (in case `x` is discontinuous at this time instant).
Since the selected step size is not close to zero, the implicit Euler step
will give a very rough approximation of the analytical integral.
A much better approximation is achieved with option
`DAE_LinearDerivativesAndConstraints` above where a step size of zero is used.

This structure is only provided for backwards compatibility.
It is numerically not reliable and should no longer be used.
"""
@enum StructureOfDAE DAE_LinearDerivativesAndConstraints DAE_ExplicitDerivatives DAE_NoSpecialStructure


"""
    simulationState = SimulationState(name, getModelResidues!, x_start,
                                      getVariableName; kwargs...)

Return a `simulationState` object that is described by a DAE
with one of the supported structures defined with enumeration [`StructureOfDAE`](@ref).

A `model` that shall be simulated with function
[`ModiaMath.simulate!`](@ref)`(model, ...)` is required to be defined as:

```julia
mutable struct ModelName <: ModiaMath.AbstractSimulationModel
    simulationState::ModiaMath.SimulationState

    # other definitions (e.g. parameters of model)
end
```

| Keyword arguments           | defaults                                           |
|:----------------------------|:---------------------------------------------------|
| structureOfDAE              | DAE_LinearDerivativesAndConstraints (see [`StructureOfDAE`](@ref)) |
| is_constraint               | fill(false, length(x_start))                       |
| has_constraintDerivatives   | false                                              |
| nz                          | 0                                                  |
| nw                          | 0                                                  |
| zDir                        | fill(0, nz)                                        |
| x_fixed                     | fill(false, length(x_start))                       |
| x_nominal                   | fill(NaN, length(x_start))                         |
| x_errorControl              | fill(true, length(x_start))                        |
| jac                         | nothing (not yet supported)                        |
| maxSparsity                 | 0.1                                                |
| hev                         | 1e-8  (only for DAE_NoSpecialStructure)            |
| scaleConstraintsAtEvents    | true  (not used; only for backwards compatibility) |
| nc                          | 0 (not used; only for backwards compatibility)     |
| getResultNames              | [`ModiaMath.getResultNames`](@ref)                 |
| storeResult!                | [`ModiaMath.storeRawResult!`](@ref)                |
| getResult                   | [`ModiaMath.getStringDictResult`](@ref)            |
| defaultTolerance            | 1e-4                                               |
| defaultStartTime            | 0.0                                                |
| defaultStopTime             | 1.0                                                |
| defaultInterval             | NaN                                                |
| defaultMaxStepSize          | NaN                                                |
| defaultMaxNumberOfSteps     | missing                                            |


# Required arguments

- `name::Union{AbstractString,Symbol}`: Name of model

- `getModelResidues!::Function`: Function with arguments `(model,t,x,derx,r,w)` to
  compute the residues `r` and auxiliary variables `w` from time `t`, vector `x` and
  its time derivative `derx`.

- `x_start::Vector{Float64}`: Start values of `x`.

- `getVariableName::Function=`[`ModiaMath.defaultVariableName`](@ref): Function that returns the
  name of a variable, given its type and its index.


# Optional (keyword) arguments:

- `structureOfDAE::`[`StructureOfDAE`](@ref): Structure of DAE.

- `is_constraint::Vector{Bool}`: = true, `residue[i]` is a constraint equation `fc` or its derivative.
                                 = false, `residue[i]` is not a constraint equation.
                                 (is only used for `structureOfDAE = DAE_LinearDerivativesAndConstraints`)

- `has_constraintDerivatives::Vector{Bool}`: if [`ModiaMath.compute_der_fc`](@ref) returns true
                        and `is_constraint[i] = true`:
                         = true , `residue[i]` is the derivative of an `fc` equation
                         = false, `residue[i]` is an `fc` equation
                        (is only used for `structureOfDAE = DAE_LinearDerivativesAndConstraints`)

- `nz::Int`: Number of event indicators

- `nw::Int`: Number of auxiliary variables (Float64 variables that are additionally computed
             and stored at communication points, and where start values can be provided
             for initialization)

- `zDir::Vector{Int}`: Interpretation of event indictors:
   zDir[i] = 0: Root is reported for both crossing directions,
           = 1: Root is reported when crossing from negative to positive direction
           = -1: Root is reported when crossing from positive to negative direction

- `x_fixed::Vector{Bool}`: = true, `x_start[i]` is fixed during initialization.
  = false, `x_start[i]` might be changed, e.g., due to an initial impulse.

- `x_nominal::Vector{Float64}`: Nominal values of `x`. if `x_nominal[i]=NaN` no nominal value
  is defined for `x[i]` and a nominal value is computed
  (`x_nominal[i] = abs(x_start[i]) > 1e-7 ? abs(x_start[i]) : 1.0`).

- `x_errorControl::Vector{Bool}`: = true, the absolute error tolerance is set to
  `0.1 * relativeTolerance * x_nominal[i]`. = false, the absolute error tolerance is
  switched off (is set to a large value). This is recommended for variables that are
  basically not limited (for example the angle of a shaft that is permantently rotating
  in the same direction and therefore becomes larger and larger).

- `hev::Float64`: Stepsize used during initialization and at event restart
  if `structureOfDAE = ModiaMath.DAE_NoSpecialStructure`.
  Otherwise `hev` is ignored.

- `scaleConstraintsAtEvents::Bool`: Dummy argument. Only kept for backwards compatibility.

- `nc::Int`: Dummy argument. Only kept for backwards compatibility.

- `jac`: Sparse Jacobian datastructure (currently not supported).

- `maxSparsity::Float64`: A sparse Jacobian is only used during simulation if
  sparseness of jac < maxSparsity (currently not supported)

- `getResultNames::Function`: Function that returns the names of the variables to be
  stored in the result data structure.

- `storeResult!::Function`: Function that stores the raw results.

- `getResult::Function`: Function that resturns the result data structure after the simulation.

- `defaultTolerance::Float64`: Model specific default relative tolerance, if not redefined in the call to
  [`ModiaMath.simulate!`](@ref).

- `defaultStartTime::Float64`: Model specific default start time in [s], if not redefined in the call to
  [`ModiaMath.simulate!`](@ref).

- `defaultStopTime::Float64`: Model specific default stop time in [s], if not redefined in the call to
  [`ModiaMath.simulate!`](@ref).

- `defaultInterval::Float64`: Model specific default interval in [s], if not redefined in the call to
  [`ModiaMath.simulate!`](@ref). Result data is stored every `defaultInterval` seconds.
  If `defaultInterval=NaN`, the default interval is computed as
  `interval = (stopTime - startTime)/500.0`.

- `defaultMaxStepSize::Float64`: Model specific default of the maximum absolute value of the step size.
  If `defaultMaxStepSize=NaN`, the default maximum step size of the integrator is used.

- `defaultMaxNumberOfSteps::Union{Int,Missing}`: Model specific default of the maximum number of
  steps to be taken by the solver in its attempt to reach the next output time.
  If `defaultMaxNumberOfSteps=missing`, the default value of the integrator is used (= 500).
"""
mutable struct SimulationState
    name::Symbol                      # Name of model
    model::Any                        # Model

    getModelResidues!::Function
    getVariableName::Function
    getResultNames::Function
    storeResult!::Function
    getResult::Function
    eventHandler::EventHandler
    structureOfDAE::StructureOfDAE
    is_constraint::Vector{Bool}
    has_constraintDerivatives::Bool

    nx::Int   # Length of x-vector
    nd::Int   # fd(der(x),x,t) = 0; nd = nx - nc
    nc::Int   #        fc(x,t) = 0; nc = count(is_constraint)
    nw::Int   # Number of auxiliary variables

    # Other model information
    nz::Int                           # Number of event indicators
    zDir::Vector{Cint}                # zDir[i] =  0: Root is reported for both crossing directions
                                     #         =  1: Root is reported when crossing from negative to positive direction
                                     #         = -1: Root is reported when crossing from positive to negative direction
    sparse::Bool                      # = true, if (sparse) jac present and shall be used
    # jac::Union{SparseMatrixCSC{Float64,Cint},Void}
                                     # Optional sparse Jacobian of DAE: der(f!,y) + cr*der(f!, yp)
                                     # (cr is a constant provided by the integrator)
    jac::Nothing
    #cg::Union{ModiaMath.SparseJacobian.ColumnGroups,Void} # if sparse, column groups of Jacobian jac
    cg::Nothing # if sparse, column groups of Jacobian jac

    # Model specific information not used by ModiaMath (can be used as default setting for ModiaMath)
    defaultTolerance::Float64    # default relative integration tolerance
    defaultStartTime::Float64    # default start time
    defaultStopTime::Float64     # default stop time
    defaultInterval::Float64     # default integration interval
    defaultMaxStepSize::Float64
    defaultMaxNumberOfSteps::Union{Int,Missing}

    # Run-time information
    time::ModiaMath.Time                        # Actual simulation time
    logger::ModiaMath.Logger                    # Logger object
    statistics::ModiaMath.SimulationStatistics  # Statistics object (complete after a full simulation run)
    tolerance::Float64                          # Relative integration tolerance
    startTime::Float64                          # Start time of the simulation
    stopTime::Float64                           # Stop time of the simulation
    interval::Float64                           # Interval of the simulation
    maxStepSize::Float64
    maxNumberOfSteps::Union{Int,Missing}

    x_start::Vector{Float64}
    x_fixed::Vector{Bool}
    x_nominal::Vector{Float64}
    x_errorControl::Vector{Bool}

    w::Vector{Float64}  # length(w) = nw

    # Auxiliary storage needed during initialization and at events
    eqInfo::ModiaMath.NonlinearEquationsInfo
    xev_old::Vector{Float64}
    xev_beforeEvent::Vector{Float64}
    xev::Vector{Float64}
    yNonlinearSolver::Vector{Float64}
    derxev::Vector{Float64}
    residues::Vector{Float64}
    residues2::Vector{Float64}
    derx_zero::Vector{Float64}
    jac2::Union{Matrix{Float64}, Nothing}
    yScale::Vector{Float64}
    rScale::Vector{Float64}
    tev::Float64
    hev::Float64
    FTOL::Float64
    scaleConstraintsAtEvents::Bool
    initialization::Bool

    # Information updated during simulation
    compute_der_fc::Bool          # = true, if derivative of constraint equations has to be returned in the residue
    storeResult::Bool             # = true, if DAE is called to store results of variables
    rawResult::ModiaMath.RawResult
    result::Any

    function SimulationState(name::Union{AbstractString,Symbol},
                             getModelResidues!::Function,
                             x_start::Vector{Float64},
                             getVariableName::Function=defaultVariableName;
                             structureOfDAE::StructureOfDAE = DAE_LinearDerivativesAndConstraints,
                             is_constraint::Vector{Bool} = fill(false, length(x_start)),
                             has_constraintDerivatives::Bool = false,
                             nc=0,
                             nz=0,
                             nw=0,
                             zDir::Vector{Int}=fill(0, nz),
                             x_fixed::Vector{Bool}=fill(false, length(x_start)),
                             x_nominal::Vector{Float64}=fill(NaN, length(x_start)),
                             x_errorControl::Vector{Bool}=fill(true, length(x_start)),
                             hev=1e-8,
                             scaleConstraintsAtEvents::Bool=true,
                             jac=nothing,
                             maxSparsity::Float64=0.1,
                             getResultNames::Function=ModiaMath.getResultNames,
                             storeResult!::Function=ModiaMath.storeRawResult!,
                             getResult::Function=ModiaMath.getStringDictResult,
                             defaultTolerance=1e-4,
                             defaultStartTime=0.0,
                             defaultStopTime=1.0,
                             defaultInterval=NaN,
                             defaultMaxStepSize=NaN,
                             defaultMaxNumberOfSteps=missing)

        # Check input arguments
        @assert(length(is_constraint) == length(x_start))
        @assert(nz >= 0)
        @assert(nw >= 0)
        @assert(length(x_fixed) == length(x_start))
        @assert(length(zDir) == nz)
        @assert(hev > 0.0)
        @assert(length(x_errorControl) == length(x_start))
        nx = length(x_start)

        if typeof(jac) != Nothing
            @assert(size(jac, 1) == nx)
            @assert(size(jac, 2) == nx)
        end

        @assert(0.0 <= maxSparsity <= 1.0)
        @assert(defaultTolerance > 0.0)
        @assert(defaultStopTime >= defaultStartTime)
        @assert(isnan(defaultInterval) || defaultInterval > 0.0)

        # Compute nominal values
        x_nominal2 = ones(nx)
        for i in 1:nx
            if isnan( x_nominal[i] )
                if abs(x_start[i]) > 1e-7
                    x_nominal2[i] = abs(x_start[i])
                end
            else
                x_nominal2[i] = x_nominal[i]
            end
        end

        yScale = ones(nx)
        for i = 1:nx
            yScale[i] = 1/x_nominal2[i]
        end

        # Compute utility elements
        nc2 = count(is_constraint)
        nd  = nx - nc2
        eventHandler = EventHandler(nz=nz)

        if structureOfDAE == DAE_LinearDerivativesAndConstraints
            eqInfo = ModiaMath.NonlinearEquationsInfo("???", nx, getEventResidues_DAE_LinearDerivativesAndConstraints!)
            jac2   = zeros(nx,nx)
            nc3    = nc2
        else
            eqInfo = ModiaMath.NonlinearEquationsInfo("???", nx, getEventResidues_DAE_NoSpecialStructure!)
            jac2   = nothing
            nc3    = missing
        end


        # If jac is provided and sparse enough, use sparse methods
        #sparse = typeof(jac) == Void ? false : nnz(jac)/(nx*nx) < maxSparsity
        #jac    = sparse ? convert(SparseMatrixCSC{Float64,Cint}, jac) : nothing
        #cg     = sparse ? ColumnGroups(jac) : nothing
        #if sparse
        #   @assert( nnz(jac) >= nx )   # otherwise jac is singular
        #end
        sparse = false
        jac = nothing
        cg = nothing

        new(Symbol(name), nothing, getModelResidues!, getVariableName, getResultNames,
            storeResult!, getResult, eventHandler, structureOfDAE, is_constraint,
            has_constraintDerivatives, nx, nd, nc2, nw, nz, zDir,
            sparse, jac, cg, defaultTolerance, defaultStartTime, defaultStopTime,
            defaultInterval, defaultMaxStepSize, defaultMaxNumberOfSteps, NaN, ModiaMath.Logger(),
            ModiaMath.SimulationStatistics(structureOfDAE, nx, nc3, sparse, sparse ? cg.ngroups : 0),
            NaN, NaN, NaN, NaN, NaN, missing,
            x_start, x_fixed, x_nominal2, x_errorControl, zeros(nw),
            eqInfo, zeros(nx), zeros(nx), zeros(nx),
            zeros(nx), zeros(nx), zeros(nx), zeros(nx), zeros(nx), jac2,
            yScale, ones(nx), 0.0, hev, 0.0,
            scaleConstraintsAtEvents, false, false, false)
    end
end


"""
    result = isNearlyEqual(x1, x2, x_nominal, tolerance)

The function returns true if x1 and x2 are nearly equal
"""
isNearlyEqual(x1::Float64, x2::Float64, x_nominal::Float64, tolerance::Float64) = abs(x1 - x2) <= max(tolerance*max(abs(x1),abs(x2)), 0.1*x_nominal*tolerance)



function getEventResidues_DAE_NoSpecialStructure!(eqInfo::ModiaMath.NonlinearEquationsInfo, y::Vector{Float64}, r::Vector{Float64})

    model = eqInfo.extraInfo
    sim::SimulationState = model.simulationState
    @assert(length(y) == eqInfo.ny)
    @assert(length(r) == eqInfo.ny)
    @assert(length(y) == sim.nx)

    residues = sim.residues
    hev      = sim.hev

    # Copy unknowns (y) to xev and derxev
    #if sim.nc == 0
    #    # Unknowns are at the left limit
    #    for i in eachindex(y)
    #        sim.xev[i]    = sim.xev_old[i]
    #        sim.derxev[i] = (sim.xev_old[i] - y[i]) / hev
    #    end
    #else # sim.nc > 0
        for i in eachindex(y)
            if sim.initialization && sim.x_fixed[i]
                # Unknowns are at the left limit
                sim.xev[i]    = sim.xev_old[i]
                sim.derxev[i] = (sim.xev_old[i] - y[i]) / hev
            else
                # Unknowns are at the right limit
                sim.xev[i]    = y[i]
                sim.derxev[i] = (y[i] - sim.xev_old[i]) / hev
            end
        end
    #end

    # Compute residues
    sim.time = sim.tev
    Base.invokelatest(sim.getModelResidues!, model, sim.tev, sim.xev, sim.derxev, residues, sim.w)

   # Copy to eq-residues and scale with h
   #=
    if sim.scaleConstraintsAtEvents
        for i = 1:sim.nd
            r[i] = residues[i]
        end

        for i = sim.nd + 1:sim.nx
            r[i] = residues[i] / hev
        end
    else
    =#
        for i = 1:sim.nx
            r[i] = residues[i]
        end
    #end
    eqInfo.lastNorm_r = norm(r, Inf)
    eqInfo.lastrScaledNorm_r = norm(sim.rScale .* r, Inf)

    # println("maxabs(r) = ", maxabs(r))
    # println("residue = "),display(r)
    # println("hev = ", hev, ", x_old = ", sim.xev_old[1], ", x = ", y[1], ", der(x) = ",sim.derxev[1], ", r = ", r[1])

    return nothing
end


function getEventResidues_DAE_LinearDerivativesAndConstraints!(eqInfo::ModiaMath.NonlinearEquationsInfo, y::Vector{Float64}, r::Vector{Float64})

    model = eqInfo.extraInfo
    sim::SimulationState = model.simulationState
    @assert(length(y) == eqInfo.ny)
    @assert(length(r) == eqInfo.ny)
    @assert(length(y) == sim.nx)

    residues  = sim.residues
    residues2 = sim.residues2

    # Copy unknowns (y) to xev and derxev
    if sim.nc == 0
        # Unknowns are at the left limit
        for i in eachindex(y)
            sim.xev[i]    = sim.xev_old[i]
            sim.derxev[i] = sim.xev_old[i] - y[i]
        end
    else # sim.nc > 0
        for i in eachindex(y)
            if sim.initialization && sim.x_fixed[i]
                # Unknowns are at the left limit
                sim.xev[i]    = sim.xev_old[i]
                sim.derxev[i] = sim.xev_old[i] - y[i]
            else
                # Unknowns are at the right limit
                sim.xev[i]    = y[i]
                sim.derxev[i] = y[i] - sim.xev_old[i]
            end
        end
    end

    # Compute residues
    sim.time = sim.tev
    Base.invokelatest(sim.getModelResidues!, model, sim.tev, sim.xev, sim.derx_zero, residues , sim.w)
    Base.invokelatest(sim.getModelResidues!, model, sim.tev, sim.xev, sim.derxev   , residues2, sim.w)

    for i = 1:sim.nx
        if sim.is_constraint[i]
            r[i] = residues[i]
        else
            r[i] = residues[i] - residues2[i]
        end
    end

    eqInfo.lastNorm_r        = norm(r, Inf)
    eqInfo.lastrScaledNorm_r = norm(sim.rScale .* r, Inf)

    # println("maxabs(r) = ", maxabs(r))
    # println("residue = "),display(r)

    return nothing
end

const small = sqrt( eps(Float64) )


function reinitialize!(model, sim::SimulationState, tev::Float64)
    nx = sim.nx
    nd = sim.nd
    nc = sim.nc

    # Determine consistent xev, der(xev)
    if ModiaMath.isLogEvents(sim)
        if sim.structureOfDAE == DAE_NoSpecialStructure
            println("        determine consistent DAE variables x,der(x) (with implicit Euler step; step size = ", sim.hev, ")")
        elseif sim.structureOfDAE == DAE_LinearDerivativesAndConstraints
            if nc == 0
                println("        for given x, determine consistent DAE variables der(x) (solving a linear equation system)")
            else
                println("        determine consistent DAE variables x,der(x) (with analytical integral over time instant)")
            end
        elseif sim.structureOfDAE == DAE_ExplicitDerivatives
            println("        for given x, compute der(x)")
        else
            error("... Cannot occur. Variable structureOfDAE = ", sim.structureOfDAE, " is wrong.")
        end
    end


    # Compute consistent xev(t+)
    if  sim.structureOfDAE == DAE_NoSpecialStructure ||
       (sim.structureOfDAE == DAE_LinearDerivativesAndConstraints && nc > 0)

        for i in eachindex(sim.xev)
            sim.xev_old[i] = sim.xev[i]
        end

        sim.tev = tev
        sim.eqInfo.lastNorm_r        = 1.0
        sim.eqInfo.lastrScaledNorm_r = 1.0

        for i in eachindex(sim.yNonlinearSolver)
            sim.yNonlinearSolver[i] = sim.xev[i]
        end

        ModiaMath.solveNonlinearEquations!(sim.eqInfo, sim.yNonlinearSolver; yScale=sim.yScale, rScale=sim.rScale, FTOL=sim.FTOL)

        # Print log message
        if ModiaMath.isLogInfos(sim)
            for i = 1:nx
                if !isNearlyEqual(sim.xev_beforeEvent[i], sim.xev[i], sim.x_nominal[i], sim.tolerance)
                    xname = sim.getVariableName(model, Category_X, i)
                    println("            ", xname, " = ", sim.xev_beforeEvent[i], " changed to ", sim.xev[i])
                end
            end
        end
    end


    # Compute consistent derivatives der(xev)(t+)
    if sim.structureOfDAE == DAE_ExplicitDerivatives
        Base.invokelatest(sim.getModelResidues!, model, sim.tev, sim.xev, sim.derx_zero, sim.derxev, sim.w)

    elseif sim.structureOfDAE == DAE_LinearDerivativesAndConstraints
        # Determine linear factors
        jac2 = sim.jac2
        sim.compute_der_fc = true
        if ModiaMath.isLogEvents(sim) && nc > 0 && sim.has_constraintDerivatives
            println("        compute der(x) with Jacobian that is constructed with model provided constraint derivatives (der(fc))")
        end
        Base.invokelatest(sim.getModelResidues!, model, sim.tev, sim.xev, sim.derx_zero, sim.residues, sim.w)

        # Determine der(fd, der(x)) and der(fc, x) using the fact that fd and der(fc) are linear in der(x)
        for j = 1:nx
            sim.derx_zero[j] = 1.0
            Base.invokelatest(sim.getModelResidues!, model, sim.tev, sim.xev, sim.derx_zero, sim.residues2, sim.w)
            for i = 1:nx
                if !sim.is_constraint[i] || sim.has_constraintDerivatives
                    # println("... 1: jac2[$i,$j]")
                    jac2[i,j] = sim.residues2[i] - sim.residues[i]
                end
            end
            sim.derx_zero[j] = 0.0
        end
        sim.compute_der_fc = false

        # Determine der(fc, x) with forward difference, if not yet computed
        if !sim.has_constraintDerivatives && nc > 0
        	inc::Float64 = 0.0
            xinc = sim.yNonlinearSolver
            for j = 1:nx
                xinc[j] = 0.0
            end
            for j = 1:nx
                xj  = sim.xev[j]
                inc = small*max( abs(xj), sim.x_nominal[j] )
                inc = (xj + inc) - xj
                sim.xev[j] += inc
                Base.invokelatest(sim.getModelResidues!, model, sim.tev, sim.xev, sim.derx_zero, sim.residues2, sim.w)
                for i = 1:nx
                    if sim.is_constraint[i]
                        # println("... 2: jac2[$i,$j]")
                        jac2[i,j] = sim.residues2[i]/inc - sim.residues[i]/inc
                    end
                end
                sim.xev[j] = xj
            end
        end

        # Solve linear system of equation: jac2*derxev = -residue
        for j = 1:nx
            sim.derxev[j] = -sim.residues[j]
        end
        ldiv!(lu!(jac2),sim.derxev)
    end
end

function eventIteration!(model, sim::SimulationState, tev::Float64)
    eh::EventHandler = sim.eventHandler

    # Initialize event iteration
    initEventIteration!(eh, tev)

    # Perform event iteration
	iter_max = 20
	success  = false
    for iter = 1:iter_max
        # Determine event branches
        for i = 1:sim.nx
            sim.xev_beforeEvent[i] = sim.xev[i]
        end
        eh.event = true
        Base.invokelatest(sim.getModelResidues!, model, tev, sim.xev, sim.derxev, sim.residues, sim.w)

        # Fix event branches and determine new consistent sim.derxev_start (and sim.xev_start during initialization)
        eh.event   = false
        eh.initial = false
        reinitialize!(model, sim, tev)

		if terminateEventIteration!(eh)
            eh.event = false
			success  = true
            break
        end
    end

	if !success
		error("Maximum number of event iterations (= $iter_max) reached")
    end
end


function initialize!(model, sim::SimulationState, t0::Float64, nt::Int, tolerance::Float64)
    eh::EventHandler = sim.eventHandler
    eh.zDir    = sim.zDir
    eh.logger  = sim.logger
    eh.initial = true
    eh.afterSimulationStart = false
    sim.eqInfo.extraInfo = model
    sim.rawResult = ModiaMath.RawResult(nt, sim.getResultNames(model))
    sim.model     = model
    sim.FTOL      = sim.tolerance    # max(sim.tolerance, eps(Float64)^(1 / 3))   # eps(Float64)^(1 / 3)   # residue tolerance for nonlinear solver
    # println("... FTOL = ", sim.FTOL)
    sim.initialization = true

    # Initialize auxiliary arrays for event iteration
    # println("... initialize!: scaleConstraintsAtEvents = ", sim.scaleConstraintsAtEvents, ", nc = ", sim.nc, ", nd = ", sim.nd)
    for i in eachindex(sim.xev),
        sim.xev[i]    = sim.x_start[i]
        sim.derxev[i] = 0.0
    end

    #if sim.scaleConstraintsAtEvents
    #    for i = 1:sim.nd
    #        sim.rScale[i] = 1.0
    #    end
    #    for i = sim.nd + 1:sim.nx
    #        sim.rScale[i] = 1.0 / sim.hev
    #    end
    #else
        for i = 1:sim.nx
            sim.rScale[i] = 1.0
        end
    #end

    # Print initial values
    if ModiaMath.isLogEvents(sim)
        println("        initial values:")
        x_table = DataFrames.DataFrame(name=String[sim.getVariableName(model, Category_X, i)  for i = 1:sim.nx],
                                       start=sim.x_start, fixed=sim.x_fixed, nominal=sim.x_nominal)

        # Print x_table, but indented and type-information removed from the heading
        io = IOBuffer()
        show(io, x_table, summary=false, rowlabel=:x)
        str=String(take!(io))
        newline = isequal('\n')
        i1=findfirst(newline, str)
        i2=findnext(newline, str, i1+1)
        i3=findnext(newline, str, i2+1)
        str2="          " * str[i1+1:i2] * str[i3+1:end]
        str3=replace(str2, "\n" => "\n          ")
        println(str3, "\n")
    end

    # Perform initial event iteration
    eventIteration!(model, sim, t0)
    sim.initialization = false

    # Compute auxiliary variables w and store all results
    computeAndStoreResult!(model, sim, t0, sim.xev, sim.derxev)
    eh.afterSimulationStart = true

    return InitInfo(sim.xev, sim.derxev;
                    y_nominal=sim.x_nominal,
                    y_errorControl=sim.x_errorControl,
                    maxTime=eh.maxTime,
                    nextEventTime=eh.nextEventTime,
                    integrateToEvent=eh.integrateToEvent,
                    terminate=eh.restart == Terminate)
end


getResidues!(model, sim::SimulationState, t::Float64, y::Vector{Float64}, yp::Vector{Float64}, residues::Vector{Float64}, hcur) =
   Base.invokelatest(sim.getModelResidues!, model, t, y, yp, residues, sim.w)


function computeAndStoreResult!(model, sim::SimulationState, t::Float64, y::Vector{Float64}, yp::Vector{Float64})
    sim.storeResult = true
    sim.time        = t
    Base.invokelatest(sim.getModelResidues!, model, t, y, yp, sim.residues, sim.w)
    sim.storeResult = false
    sim.storeResult!(sim.rawResult, model, t, y, yp, sim.w)
    return nothing
end


function terminate!(model, sim::SimulationState, t::Float64, y::Vector{Float64}, yp::Vector{Float64})
    eh::EventHandler = sim.eventHandler
    eh.terminal = true
    sim.time    = t
    Base.invokelatest(sim.getModelResidues!, model, t, y, yp, sim.residues, sim.w)
    eh.terminal = false
    #sim.result,nt = Result.getDictResult(sim.rawResult)
    sim.result, nt  = sim.getResult(model, sim.rawResult)
    ModiaMath.set_nResultsForSimulationStatistics!(sim.statistics, nt)
    return sim.result
end


function processEvent!(model, sim::SimulationState, t::Float64, y::Vector{Float64}, yp::Vector{Float64}, eventInfo::EventInfo)
    eh::EventHandler = sim.eventHandler

    # Event iteration
    for i in eachindex(sim.xev)
        sim.xev[i]    = y[i]
        sim.derxev[i] = yp[i]
    end
    eventIteration!(model, sim, t)
    for i in eachindex(sim.xev)
        y[i]  = sim.xev[i]
        yp[i] = sim.derxev[i]
    end

    # Compute auxiliary variables w and store all results
    computeAndStoreResult!(model, sim, t, y, yp)

    eventInfo.restart          = eh.restart
    eventInfo.maxTime          = eh.maxTime
    eventInfo.nextEventTime    = eh.nextEventTime
    eventInfo.integrateToEvent = eh.integrateToEvent

    return nothing
end


function getEventIndicators!(model, sim::SimulationState, t::Float64, y::Vector{Float64}, yp::Vector{Float64}, z::Vector{Float64})
    eh::EventHandler = sim.eventHandler
    eh.crossing = true
    sim.time    = t
    Base.invokelatest(sim.getModelResidues!, model, t, y, yp, sim.residues, sim.w)
    eh.crossing = false

    for i = 1:eh.nz
        z[i] = eh.z[i]
    end
    return nothing
end



#=
"""
    computeFullJacobian!(model, sim, t, y, yp, residues, cj, ewt)

Compute full Jacobian

```julia
jac = der(f,y) + cj*der(f,yp).
```

numerically by finite differences
"""
function computeJacobian!(model, sim::SimulationState, t::Float64, y::Vector{Float64}, yp::Vector{Float64}, r::Vector{Float64},
                          fulljac::Matrix{Float64}, hcur::Float64, cj::Float64, ewt::Vector{Float64})
    # Compute one column of the Jacobian once at a time
    inc::Float64 = 0.0
    yj::Float64  = 0.0
    ypj::Float64 = 0.0
    for j = 1:sim.nx
        # Determine increment for column j
        yj  = y[j]
        ypj = yp[j]
        inc = max( small*max( abs(yj), abs(hcur*ypj) ), 1.0/ewt[j] )
        if hcur*ypj < 0.0
            inc = -inc
        end
        inc = (yj + inc) - yj
        y[j]  += inc
        yp[j] += cj*inc

        # Compute residues
        Base.invokelatest(sim.getModelResidues!, model, t, y, yp, sim.residues, sim.w)

        # Compute finite differences and store them in the j column of the Jacobian
        for i = 1:sim.nx
            fulljac[i,j] = sim.residues[i]/inc - r[j]/inc
        end

        # Reset y and yp
        y[j]  = yj
        yp[j] = ypj
    end
end
=#
