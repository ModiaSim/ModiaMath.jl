# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.DAE (ModiaMath/DAE/_module.jl)
#

# Utility functions used below
defaultVariableName(model::Any, vcat::VariableCategory, vindex::Int) = vcat == Category_X ? "x[" * string(vindex) * "]" :
                                                                      (vcat == Category_W ? "w[" * string(vindex) * "]" :
                                                                                        "der(x[" * string(vindex) * "])" )

#nameVector(name, nNames::Int) = [Symbol(name, "[", i, "]") for i = 1:nNames]                                                                            
#fcNameVector(fc, names::AbstractVector) = [ Symbol(fc, "(", names[i], ")") for i in eachindex(names)]

nameVector(name, nNames::Int) = [string(name, "[", i, "]") for i = 1:nNames]                                                                            
fcNameVector(fc, names::AbstractVector) = [ string(fc, "(", names[i], ")") for i in eachindex(names)]

                                                                    
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


function getResultNames(model::Any) 
    sim::SimulationState = model.simulationState
    resultNames = ["time"; String[sim.getVariableName(model, Category_X, i)  for i = 1:sim.nx];
                          String[sim.getVariableName(sim, Category_DERX, i) for i = 1:sim.nx];
                          String[sim.getVariableName(sim, Category_W, i)    for i = 1:sim.nw] ]
    return resultNames
end


"""
    simulationState = SimulationState(
          name, getModelResidues!, x_start, 
          getVariableName=ModiaMath.defaultVariableName;
          nc=0, nz=0, nw=0,
          zDir                     = fill(0, nz),                   
          x_fixed                  = fill(false, length(x_start)),   
          x_nominal                = fill(NaN, length(x_start)),   
          x_errorControl           = fill(true, length(x_start)),
          w_start                  = fill(NaN,nw),
          w_fixed                  = fill(false,nw),
          w_nominal                = fill(NaN, length(w_start)),   
          jac                      = nothing, 
          maxSparsity              = 0.1,
          hev                      = 1e-8,
          scaleConstraintsAtEvents = true,               
          getResultNames::Function = ModiaMath.getResultNames, 
          storeResult!::Function   = ModiaMath.storeRawResult!,
          getResult::Function      = ModiaMath.getStringDictResult,
          defaultTolerance         = 1e-4, 
          defaultStartTime         = 0.0,
          defaultStopTime          = 1.0, 
          defaultInterval          = NaN)
                     
Return a `simulationState` object. A `model` that shall be simulated with function
[`ModiaMath.simulate!`](@ref)`(model, ...)` is required to be defined as:

```julia
mutable struct ModelName <: ModiaMath.AbstractSimulationModel
    simulationState::ModiaMath.SimulationState

    # other definitions (e.g. parameters of model)
end
```


# Required arguments

- `name::Union{AbstractString,Symbol}`: Name of model

- `getModelResidues!::Function`: Function with arguments `(model,t,x,derx,r,w)` to
  compute the residues `r` and auxiliary variables `w` from time `t`, vector `x` and
  its time derivative `derx`.
 
- `x_start::Vector{Float64}`: Start values of `x`.

- `getVariableName::Function=ModiaMath.defaultVariableName`: Function that returns the
  name of a variable, given its type and its index.


# Optional (keyword) arguments:

- `nc::Int`: Number of constraints functions (= length(fc))

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

- `w_start::Vector{Float64}`: Start values for auxiliary variables `w`.
  If `w_start[i] = NaN`, then no start value for `w[i]` is defined and
  `w[i]` is ignored during initialization.
  If `w_start[i] != NaN`, an initial equation `w[i] = w_start[i]` is utilized
  during initialization.

- `w_fixed::Vector{Bool}`: = true (and `w_start[i] != NaN`), `w_start[i]` is fixed 
  during initialization. = false (and `w_start[i] != NaN`, `w_start[i]`
  might be changed, e.g., due to an initial impulse.

- `w_nominal::Vector{Float64}`: Nominal values of `w`. if `w_nominal[i]=NaN` no nominal value
  is defined for `w[i]` and a nominal value is computed
  (`w_nominal[i] = abs(w_start[i]) > 1e-7 ? abs(w_start[i]) : 1.0`).

- `hev::Float64`: Stepsize used during initialization and at event restart.

- `scaleConstraintsAtEvents::Bool`: = true, constraint equations are scaled during 
  initialization and at event restart (currently, this setting is ignored).

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

    nx::Int   # Length of x-vector
    nd::Int   # fd(der(x),x,t) = 0; nd = length(fd)
    nc::Int   #        fc(x,t) = 0; nc = length(fc)
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

    # Run-time information 
    time::ModiaMath.Time                        # Actual simulation time
    logger::ModiaMath.Logger                    # Logger object
    statistics::ModiaMath.SimulationStatistics  # Statistics object (complete after a full simulation run)
    tolerance::Float64                          # Relative integration tolerance
    startTime::Float64                          # Start time of the simulation
    stopTime::Float64                           # Stop time of the simulation
    interval::Float64                           # Interval of the simulation

    x_start::Vector{Float64}
    x_fixed::Vector{Bool}   
    x_nominal::Vector{Float64}
    x_errorControl::Vector{Bool}     

    w::Vector{Float64}  # length(w) = nw
    w_start::Vector{Float64}
    w_fixed::Vector{Bool}
    w_nominal::Vector{Float64}
    w_start_used::Vector{Int}                # The indices of w_start, where w_start[i] != NaN.
                                             # Only these w_start values are used for initialization
        
    # Auxiliary storage needed during initialization and at events
    eqInfo::ModiaMath.NonlinearEquationsInfo
    xev_old::Vector{Float64}
    xev_beforeEvent::Vector{Float64}
    xev::Vector{Float64}
    yNonlinearSolver::Vector{Float64}
    derxev::Vector{Float64}
    residues::Vector{Float64}
    yScale::Vector{Float64}
    rScale::Vector{Float64}
    tev::Float64
    hev::Float64
    FTOL::Float64
    scaleConstraintsAtEvents::Bool
    initialization::Bool
      
    # Information updated during simulation   
    storeResult::Bool             # = true, if DAE is called to store results of variables
    rawResult::ModiaMath.RawResult
    result::Any
   
    function SimulationState(name::Union{AbstractString,Symbol},
                            getModelResidues!::Function, 
                            x_start::Vector{Float64},
                            getVariableName::Function=defaultVariableName;
                            nc=0,
                            nz=0,
                            nw=0,
                            zDir::Vector{Int}=fill(0, nz),                   
                            x_fixed::Vector{Bool}=fill(false, length(x_start)),   
                            x_nominal::Vector{Float64}=fill(NaN, length(x_start)),                                                   
                            x_errorControl::Vector{Bool}=fill(true, length(x_start)),
                            w_start::Vector{Float64}=fill(NaN,nw),
                            w_fixed::Vector{Bool}=fill(false,nw),
                            w_nominal::Vector{Float64}=fill(NaN, length(w_start)),
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
                            defaultInterval=NaN)
  
        # Check input arguments                                        
        @assert(nc >= 0)
        @assert(nc <= length(x_start))
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

        w_nominal2 = ones(nw)
        for i in 1:nw
            if isnan( w_nominal[i] )
                if abs(w_start[i]) > 1e-7
                    w_nominal2[i] = abs(w_start[i])
                end
            else
                w_nominal2[i] = w_nominal[i]
            end
        end

        # Compute w_start values that are used during initialization
        w_start_used = Float64[]
        for i in eachindex(w_start)
            if !isnan(w_start[i])
                push!(w_start_used, i)
            end
        end

        # Compute utility elements
        nd = nx - nc
        eventHandler = EventHandler(nz=nz)
            
        eqInfo = ModiaMath.NonlinearEquationsInfo("???", nx, getEventResidues!)

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
            storeResult!, getResult, eventHandler, nx, nd, nc, nw, nz, zDir,
            sparse, jac, cg, defaultTolerance, defaultStartTime, defaultStopTime,
            defaultInterval, NaN, ModiaMath.Logger(),
            ModiaMath.SimulationStatistics(nx, sparse, sparse ? cg.ngroups : 0),
            NaN, NaN, NaN, NaN,
            x_start, x_fixed, x_nominal2, x_errorControl, zeros(nw),  
            w_start, w_fixed, w_nominal2, w_start_used,
            eqInfo, zeros(nx), zeros(nx), zeros(nx),
            zeros(nx), zeros(nx), zeros(nx), yScale, ones(nx), 0.0, hev, 0.0, 
            scaleConstraintsAtEvents, false, false)
    end
end


"""
    result = isNearlyEqual(x1, x2, x_nominal, tolerance)
    
The function returns true if x1 and x2 are nearly equal
"""
isNearlyEqual(x1::Float64, x2::Float64, x_nominal::Float64, tolerance::Float64) = abs(x1 - x2) <= max(tolerance*max(abs(x1),abs(x2)), 0.1*x_nominal*tolerance)



function getEventResidues!(eqInfo::ModiaMath.NonlinearEquationsInfo, y::Vector{Float64}, r::Vector{Float64})
    #=
       DAE:
           w = fw(der(x),x,t)             
           0 = fd(der(x),x,t,w)          [der(fd,der(x));
           0 = fc(x,t,w)                  der(fc,x)]      is regular (under the assumption that w is inlined and
                                                                      w[i] used in fc(...) do not depend on der(x))

       Initialization/re-initialization equation:
           0 = f(y)  ->  solve for y

       Re-Initialization at event:
           # It is assumed that x might have been changed at the event instant.
           # For example when an impulse occured: v(after_event) = -eps*v(before_event).
           # However, it is implicitely assumed (this is not checked) that a Dirac impulse
           # only occurs in der(x) and not in x. This in turn means that if ni is the highest
           # occuring derivative of x[i], then at most the n-1-derivative of x[i] is allowed to
           # be discontinuous and all lower derivatives must be continuous at the event. 
           #     If this requirement is fulfilled, it can be guaranteed that integrating over
           # the discontinuity will yield a mathematically well-defined solution.
           # If this requirement is not fulfilled, it is unclear whether a mathematically
           # well-defined solution exists and most likely a somewhat arbitrary value is computed.
           
           dim(fc) = 0:
               # Any "x" can be used (note, it is assumed that der(fd,der(x)) is regular).
               # Therefore, the right limit of x is fixed and derivatives der(x) could be used
               # as iteration variables. However, since there is no nominal value for der(x),
               # instead the following approximation is used where the iteration variables y[i]
               # are the virtual previous values of x
                        x[i] = x_ev[i]                     # x_ev the value of x after the event
                   der(x[i]) = (x_ev[i] - y[i]) / hev

               # If der(x[i]) can be explicitely computed, then there is no iteration variable
               # associated with der(x[i]) or the virtual previous value of x. So vector y has one 
               # element less (If all der(x[i]) can be explicitly computed (so ODE), then no equation
               # must be solved at all):
                   x[i] = x_ev[i]
                   der(x[i]) is explicitely computed
 
               # If fd is linear in der(x), which seems to be always the case for physical systems,
               # no nonlinear equation must be solved. Instead, the Jacobian with respect to der(x) 
               # is computed and then a linear equation is solved to compute der(x)
               # (this is not yet done):
                   0 = Jd(x,t)*der(x) + hd(x,t)  -> der(x) = solve(Jd, -hd)
               
                   # Determination of Jd and hd:
                       -      hd = fd(der(x)=0  , x, t)
                       - Jd[:,i] = fd(der(x)=e_i, x, t) - hd(x,t)  

               # If fd is linear in der(x) and some der(x[i]) can be explicitely computed, then the
               # linear equation system to be solved for is smaller. Furthermore, the Jacobian could
               # also be sparse. There are several options, such as:
               #   (a) Only support sparse der(x)-Jacobians. In this case explicitely solvable
               #       derivatives are included as a special case. The model must provide the structure 
               #       of the Jacobian.
               #   (b) The generated code solves already explicitely for the derivatives at an event,
               #       so from the outside, all derivatives are explicitely computed.
           
           dim(fc) > 0:
               # fc(x,t) might be non-zero at an event (for example, if an impulse occured
               # and v(after_event) = -eps*v(before_event), then the velocity constraint of a 
               # closed loop system might be no longer fulfilled.
               # Therefore, the left limit of x is fixed and the right limit are the iteration variables y
               # (that are determined, so that fc(x,t) = 0 and fd(der(x),x,t)=0)
                        x[i] = y[i]                     # the right limit of x
                   der(x[i]) = (y[i] - x_ev[i]) / hev   # the derivative at the right limit of x

               # If der(x[i]) can be explicitely computed, then x[i] might still need 
               # to be modified, such that fc(x,t) = 0. For this reason, the dimension of the 
               # reinitialization problem is not changed, just der(x[i]) is analytically computed
               # and the residue is computed as the difference to the numerically computed derivative.
               # Therefore, no essential simplification occurs when utilizing this feature, just
               # the residue computation is more robust.

               # If fd is linear in der(x), still the nonlinear equation fc(x,t)=0 must be solved.
               # Again no simplification for the re-initialization algorithm seems to be possible,
               # if it is known that fd is linear in the derivative.


        Initialization at t_start:
           # Similar to re-initialization, but additionally with x_fixed[i] = true it is defined
           # that x[i] is fixed and should not be changed.
           # Furthermore, w_fixed[i] 
           if x_fixed[i] or dim(fc) = 0
              # Right limit of x is fixed, left limit of x are the iteration variables y
                   x[i] = x_start[i]
              der(x[i]) = (x_start[i] - y[i]) / hev

              # If dim(fc) = 0, a linear equation can be solved to compute der(x)
           else
              # Left limit of x is fixed, right limit of x are the iteration variables y
                   x[i] = y[i]
              der(x[i]) = (y[i] - x_start[i]) / hev
           end
    =#

    model = eqInfo.extraInfo
    sim::SimulationState = model.simulationState
    @assert(length(y) == eqInfo.ny)
    @assert(length(r) == eqInfo.ny)
    @assert(length(y) == sim.nx)
   
    residues = sim.residues
    hev      = sim.hev
   
    # Copy unknowns (y) to xev and derxev
    if sim.nc == 0
        # Unknowns are at the left limit
        for i in eachindex(y)
            sim.xev[i]    = sim.xev_old[i]
            sim.derxev[i] = (sim.xev_old[i] - y[i]) / hev            
        end 
    else # sim.nc > 0
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
    end
   
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


function reinitialize!(model, sim::SimulationState, tev::Float64)
    nx = sim.nx   
    nd = sim.nd
    nc = sim.nc
 
    # Determine consistent initial conditions
    for i in eachindex(sim.xev)
        sim.xev_old[i] = sim.xev[i]   
    end
  
    # Determine consistent xev, der(xev)
    if ModiaMath.isLogEvents(sim)
        if nc == 0
            println("        for given x, determine consistent DAE variables der(x) (with implicit Euler step; step size = ", sim.hev, ")")
        else
            println("        determine consistent DAE variables x,der(x) (with implicit Euler step; step size = ", sim.hev, ")")
        end
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


function eventIteration!(model, sim::SimulationState, tev::Float64)
    eh::EventHandler = sim.eventHandler 

    # Initialize event iteration
    initEventIteration!(eh, tev)

    # Perform event iteration
    while true
        # Determine event branches
        for i = 1:sim.nx
            sim.xev_beforeEvent[i] = sim.xev[i]
        end
        eh.event = true
        Base.invokelatest(sim.getModelResidues!, model, tev, sim.xev, sim.derxev, sim.residues, sim.w)
      
        if terminateEventIteration!(eh)
            eh.event = false
            break
        end
      
        # Fix event branches and determine new consistent sim.derxev_start (and sim.xev_start during initialization)
        eh.event   = false
        eh.initial = false
        reinitialize!(model, sim, tev)
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
    sim.FTOL      = max(sim.tolerance, eps(Float64)^(1 / 3))   # eps(Float64)^(1 / 3)   # residue tolerance for nonlinear solver
    # println("... FTOL = ", sim.FTOL)
    sim.initialization = true

    # Initialize auxiliary arrays for event iteration   
    # println("... initialize!: scaleConstraintsAtEvents = ", sim.scaleConstraintsAtEvents, ", nc = ", sim.nc, ", nd = ", sim.nd)
    for i in eachindex(sim.xev),
        sim.xev[i]    = sim.x_start[i]
        sim.derxev[i] = 0.0
    end

    if sim.scaleConstraintsAtEvents
        for i = 1:sim.nd
            sim.rScale[i] = 1.0
        end
        for i = sim.nd + 1:sim.nx
            sim.rScale[i] = 1.0 / sim.hev
        end
    else
        for i = 1:sim.nx
            sim.rScale[i] = 1.0
        end
    end

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




const small = sqrt( eps(Float64) )


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

