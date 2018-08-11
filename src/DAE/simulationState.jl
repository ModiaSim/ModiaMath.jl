# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.DAE (ModiaMath/DAE/_module.jl)
#

# Utility functions used below
defaultVariableName(model::Any, vcat::VariableCategory, vindex::Int) = vcat == Category_X ? "x["*string(vindex)*"]" :
                                                                      (vcat == Category_W ? "w["*string(vindex)*"]" :
                                                                                        "der(x["*string(vindex)*"])" )

#nameVector(name, nNames::Int) = [Symbol(name, "[", i, "]") for i = 1:nNames]                                                                            
#fcNameVector(fc, names::AbstractVector) = [ Symbol(fc, "(", names[i], ")") for i in eachindex(names)]

nameVector(name, nNames::Int) = [string(name, "[", i, "]") for i = 1:nNames]                                                                            
fcNameVector(fc, names::AbstractVector) = [ string(fc, "(", names[i], ")") for i in eachindex(names)]

                                                                    
function getVariableName(model::Any,vcat::VariableCategory, vindex::Int, nx=0;
                         xNames::Vector{String}    = nameVector("x",nx), 
                         derxNames::Vector{String} = fcNameVector("der",xNames),
                         wNames::Vector{String}    = String[])
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
   resultNames = ["time"; String[sim.getVariableName(model,Category_X,i)  for i=1:sim.nx];
                          String[sim.getVariableName(sim,Category_DERX,i) for i=1:sim.nx];
                          String[sim.getVariableName(sim,Category_W,i)    for i=1:sim.nw] ]
   return resultNames
end



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
   #jac::Union{SparseMatrixCSC{Float64,Cint},Void}
                                     # Optional sparse Jacobian of DAE: der(f!,y) + cr*der(f!, yp)
                                     # (cr is a constant provided by the integrator)
   @static if VERSION >= v"0.7.0-DEV.2005"
       jac::Nothing
       #cg::Union{ModiaMath.SparseJacobian.ColumnGroups,Void} # if sparse, column groups of Jacobian jac
       cg::Nothing # if sparse, column groups of Jacobian jac
   else
       jac::Void
       #cg::Union{ModiaMath.SparseJacobian.ColumnGroups,Void} # if sparse, column groups of Jacobian jac
       cg::Void # if sparse, column groups of Jacobian jac

   end
   # Model specific information not used by ModiaMath (can be used as default setting for ModiaMath)
   defaultTolerance::Float64    # default relative integration tolerance
   defaultStartTime::Float64    # default start time
   defaultStopTime::Float64     # default stop time
   defaultInterval::Float64     # default integration interval

   # Run-time information 
   time::ModiaMath.Time                        # Actual simulation time
   logger::ModiaMath.Logger                    # Logger object
   statistics::ModiaMath.SimulationStatistics  # Statistics object (complete after a full simulation run)
   tolerance::Float64                          # Telative integration tolerance
   startTime::Float64                          # Start time of the simulation
   stopTime::Float64                           # Stop time of the simulation
   interval::Float64                           # Interval of the simulation

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
   rScale::Vector{Float64}
   tev::Float64
   hev::Float64
   FTOL::Float64
   scaleConstraintsAtEvents::Bool
   initialization::Bool
   use_x_fixed::Bool
      
   # Information updated during simulation   
   storeResult::Bool             # = true, if DAE is called to store results of variables
   rawResult::ModiaMath.RawResult
   result::Any
   
   function SimulationState(name,
                            getModelResidues!::Function, 
                            x_start::Vector{Float64},
                            getVariableName::Function = defaultVariableName;
                            nc = 0,
                            nz = 0,
                            nw = 0,
                            zDir::Vector{Int} = fill(0,nz),                   
                            x_fixed::Vector{Bool}        = fill(false,length(x_start)),                                                      
                            x_errorControl::Vector{Bool} = fill(true,length(x_start)),
                            hev = 1e-8,
                            scaleConstraintsAtEvents::Bool = true,               
                            jac=nothing, 
                            maxSparsity::Float64=0.1,
                            getResultNames::Function = ModiaMath.getResultNames, 
                            storeResult!::Function = ModiaMath.storeRawResult!,
                            getResult::Function = ModiaMath.getStringDictResult,
                            defaultTolerance=1e-4, 
                            defaultStartTime=0.0,
                            defaultStopTime=1.0, 
                            defaultInterval=(defaultStopTime-defaultStartTime)/500.0)
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
      @static if VERSION >= v"0.7.0-DEV.2005"
       if typeof(jac) != Nothing
         @assert(size(jac,1) == nx)
         @assert(size(jac,2) == nx)
       end
      else
       if typeof(jac) != Void
         @assert(size(jac,1) == nx)
         @assert(size(jac,2) == nx)
       end
      end
      @assert(0.0 <= maxSparsity <= 1.0)
      @assert(defaultTolerance > 0.0)
      @assert(defaultStopTime >= defaultStartTime)
      @assert(defaultInterval > 0.0)

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
          x_start, x_fixed, ones(nx), x_errorControl, zeros(nw),  
          eqInfo, zeros(nx), zeros(nx), zeros(nx),
          zeros(nx), zeros(nx), zeros(nx), ones(nx), 0.0, hev, 0.0, 
          scaleConstraintsAtEvents, false, false, false)
   end
end


"""
    result = isNearlyEqual(x1,x2)
    
The function returns true if x1 and x2 are nearly equal
"""
isNearlyEqual(x1::Float64, x2::Float64) = abs(x1) > 1e-8 ? abs(x1-x2)/max(abs(x1)) <= 1e-8 : abs(x1-x2) <= 1e-8


function getEventResidues!(eqInfo::ModiaMath.NonlinearEquationsInfo, y::Vector{Float64}, r::Vector{Float64})
   model = eqInfo.extraInfo
   sim::SimulationState = model.simulationState
   @assert(length(y) == eqInfo.ny)
   @assert(length(r) == eqInfo.ny)
   @assert(length(y) == sim.nx)
   
   residues = sim.residues
   hev      = sim.hev
   
   # Copy unknowns (y) to xev and derxev
   for i in eachindex(y)
      if sim.use_x_fixed && sim.x_fixed[i]
         # Unknowns are at the left limit
         sim.xev[i]    = sim.xev_old[i]
         sim.derxev[i] = (sim.xev_old[i] - y[i])/hev           
      else
         # Unknowns are at the right limit
         sim.xev[i]    = y[i]
         sim.derxev[i] = (y[i] - sim.xev_old[i])/hev
      end
   end
   
   # Compute residues
   sim.time = sim.tev
   Base.invokelatest(sim.getModelResidues!, model, sim.tev, sim.xev, sim.derxev, residues, sim.w)  
                              
   # Copy to eq-residues and scale with h
   if sim.scaleConstraintsAtEvents
      for i = 1:sim.nd
         r[i] = residues[i]
      end
      for i = sim.nd+1:sim.nx
         r[i] = residues[i]/hev
      end
   else
      for i = 1:sim.nx
         r[i] = residues[i]
      end
   end
   eqInfo.lastNorm_r = norm(r,Inf)
   eqInfo.lastrScaledNorm_r = norm(sim.rScale.*r,Inf)
   
   # println("maxabs(r) = ", maxabs(r))
   # println("residue = "),display(r)
   # println("hev = ", hev, ", x_old = ", sim.xev_old[1], ", x = ", y[1], ", der(x) = ",sim.derxev[1], ", r = ", r[1])
   
   return nothing
end


function reinitialize!(model, sim::SimulationState, tev::Float64, xevIsConsistent::Bool)
   nx = sim.nx   
   nd = sim.nd
   nc = sim.nc
 
   # Determine consistent initial conditions
   for i in eachindex(sim.xev)
      sim.xev_old[i] = sim.xev[i]   
   end
   if !xevIsConsistent
      # xev is potentially not consistent; determine consistent xev
      if ModiaMath.isLogEvents(sim)
         println("        determine consistent DAE variables x (with implicit Euler step; step size = ",sim.hev,")")
         if sim.initialization
            println("            (fixed attribute defines whether unknowns are left or right limit values)")
         end
      end
      sim.tev = tev
      sim.eqInfo.lastNorm_r        = 1.0
      sim.eqInfo.lastrScaledNorm_r = 1.0      
      for i in eachindex(sim.yNonlinearSolver)
         sim.yNonlinearSolver[i] = sim.xev[i]
      end
      if sim.initialization
         sim.use_x_fixed = true
      end

      ModiaMath.solveNonlinearEquations!(sim.eqInfo, sim.yNonlinearSolver; rScale = sim.rScale, FTOL = sim.FTOL)
      sim.use_x_fixed = false
      
      # Copy result to sim.xev
      for i = 1:nx
         sim.xev[i] = sim.yNonlinearSolver[i]
      end
      
      # Only during initialization: Copy fixed = true values
      if sim.initialization
         for i = 1:nx
            if sim.x_fixed[i] 
               sim.xev[i] = sim.xev_old[i]
            end
         end
      end
      
      # Copy xev to xev_old    
      for i = 1:nx
         sim.xev_old[i] = sim.xev[i]
      end
   end

   # Print log message
   if ModiaMath.isLogInfos(sim)
      for i = 1:nx
         if !isNearlyEqual(sim.xev_beforeEvent[i], sim.xev[i])
            xname = sim.getVariableName(model,Category_X, i)
            println("            ", xname, " = ", sim.xev_beforeEvent[i], " changed to ", sim.xev[i])
         end
      end
   end
   
   # Determine consistent derxev
   sim.tev = tev
   sim.eqInfo.lastNorm_r        = 1.0
   sim.eqInfo.lastrScaledNorm_r = 1.0  
   for i in eachindex(sim.yNonlinearSolver)
      sim.yNonlinearSolver[i] = sim.xev[i]
   end   
   sim.use_x_fixed = false
   if ModiaMath.isLogEvents(sim)
      println("        determine consistent DAE variables der(x) (with implicit Euler step; step size = ", sim.hev, ")")
   end
      
   ModiaMath.solveNonlinearEquations!(sim.eqInfo, sim.yNonlinearSolver; rScale = sim.rScale, FTOL = sim.FTOL) 
end


function eventIteration!(model, sim::SimulationState, tev::Float64)
   eh::EventHandler = sim.eventHandler 
   xevIsConsistent::Bool   = false

   # Initialize event iteration
   initEventIteration!(eh,tev)

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
      if eh.initial
         eh.initial = false
         xevIsConsistent = sim.nc == 0
         if sim.nc > 0
            if norm(sim.residues[sim.nd+1:end], Inf) < 1e-3
               xevIsConsistent = true
            end
         end
      else
         xevIsConsistent = true
         for i = 1:sim.nx
            if sim.xev[i] != sim.xev_beforeEvent[i]
               xevIsConsistent = false
               break
            end
         end
      end
      
      # Fix event branches and determine new consistent sim.xev_start, sim.derxev_start
      eh.event = false
      reinitialize!(model, sim, tev, xevIsConsistent)
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
   sim.FTOL      =  eps(Float64)^(1/3)   # residue tolerance for nonlinear solver
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
      for i = sim.nd+1:sim.nx
         sim.rScale[i] = 1.0/sim.hev
      end
   else
      for i = 1:sim.nx
         sim.rScale[i] = 1.0
      end
   end
   # println("... initialize!: rScale = ", sim.rScale)
   
   # Perform initial event iteration
   eventIteration!(model, sim, t0)
   sim.initialization = false   

   # Compute auxiliary variables w and store all results
   computeAndStoreResult!(model, sim, t0, sim.xev, sim.derxev)
   eh.afterSimulationStart = true

   return InitInfo(sim.xev, sim.derxev;
                   y_nominal        = sim.x_nominal,
                   y_errorControl   = sim.x_errorControl,
                   maxTime          = eh.maxTime,
                   nextEventTime    = eh.nextEventTime,
                   integrateToEvent = eh.integrateToEvent,
                   terminate        = eh.restart==Terminate)
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
   sim.result,nt  = sim.getResult(model, sim.rawResult)
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
