# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.DAE (ModiaMath/DAE/_module.jl)
#

# Functions that can be called in getResidues! function

""" 
    ModiaMath.getTime(m::ModiaMath.[AbstractSimulationModel|SimulationState])

Return actual **simulation time**.
"""
getTime(sim::SimulationState)                 = sim.time
getTime(m::ModiaMath.AbstractSimulationModel) = m.simulationState.time



""" 
    ModiaMath.getStartTime(m::ModiaMath.[AbstractSimulationModel|SimulationState])

Return **startTime** of the actual simulation run.
"""
getStartTime(sim::SimulationState)                 = sim.startTime
getStartTime(m::ModiaMath.AbstractSimulationModel) = m.simulationState.startTime


""" 
    ModiaMath.getStopTime(m::ModiaMath.[AbstractSimulationModel|SimulationState])

Return **stopTime** of the actual simulation run (when the simulation will be terminated).
"""
getStopTime(sim::SimulationState)                 = sim.stopTime
getStopTime(m::ModiaMath.AbstractSimulationModel) = m.simulationState.stopTime


""" 
    ModiaMath.getTolerance(m::ModiaMath.[AbstractSimulationModel|SimulationState])

Return (relative) **tolerance** of the actual simulation run.
"""
getTolerance(sim::SimulationState)                 = sim.tolerance
getTolerance(m::ModiaMath.AbstractSimulationModel) = m.simulationState.tolerance


""" 
    ModiaMath.get_is_constraint(m::ModiaMath.[AbstractSimulationModel|SimulationState])

Return reference to `is_constraint` vector (this vector can be modified in `getModelResidues!` if
[`ModiaMath.isEvent`](@ref) returns true.
"""
get_is_constraint(sim::SimulationState)                 = sim.is_constraint
get_is_constraint(m::ModiaMath.AbstractSimulationModel) = m.simulationState.is_constraint


""" 
    ModiaMath.compute_der_fc(m::ModiaMath.[AbstractSimulationModel|SimulationState])

If true is returned, return the derivative of the constraint equation in residues `r[i]`. 
Otherwise return constraint equations.
"""
compute_der_fc(sim::SimulationState)                 = sim.compute_der_fc
compute_der_fc(m::ModiaMath.AbstractSimulationModel) = m.simulationState.compute_der_fc


"""
    ModiaMath.isInitial(m::ModiaMath.[AbstractSimulationModel|SimulationState])

Return true, if **initialization phase** of simulation.
"""
isInitial(sim::SimulationState)                 = sim.eventHandler.initial
isInitial(m::ModiaMath.AbstractSimulationModel) = m.simulationState.eventHandler.initial



"""
    ModiaMath.isTerminal(m::ModiaMath.[AbstractSimulationModel|SimulationState])

Return true, if **terminal phase** of simulation.
"""
isTerminal(sim::SimulationState)                 = sim.eventHandler.terminal
isTerminal(m::ModiaMath.AbstractSimulationModel) = m.simulationState.eventHandler.terminal



"""
    ModiaMath.isEvent(m::ModiaMath.[AbstractSimulationModel|SimulationState])

Return true, if **event phase** of simulation (including initialization)
"""
isEvent(sim::SimulationState)                 = sim.eventHandler.event
isEvent(m::ModiaMath.AbstractSimulationModel) = m.simulationState.eventHandler.event



"""
    ModiaMath.isZeroCrossing(m::ModiaMath.[AbstractSimulationModel|SimulationState])

Return true, if simulation model shall compute zero-crossing functions (as required by the integrator).
"""
isZeroCrossing(sim::SimulationState)                 = sim.eventHandler.crossing
isZeroCrossing(m::ModiaMath.AbstractSimulationModel) = m.simulationState.eventHandler.crossing



"""
    ModiaMath.isAfterSimulationStart(m::ModiaMath.[AbstractSimulationModel|SimulationState])

Return true, if after start of simulation and false if during initialization.
"""
isAfterSimulationStart(sim::SimulationState)                 = sim.eventHandler.afterSimulationStart
isAfterSimulationStart(m::ModiaMath.AbstractSimulationModel) = m.simulationState.eventHandler.afterSimulationStart



"""
    ModiaMath.isStoreResult(m::ModiaMath.[AbstractSimulationModel|SimulationState])

Return true, if **communication point** and variables shall be stored in the result data structure.
"""
isStoreResult(sim::SimulationState)                 = sim.storeResult
isStoreResult(m::ModiaMath.AbstractSimulationModel) = m.simulationState.storeResult


ModiaMath.Logging.isLogInfos(sim::SimulationState) = ModiaMath.Logging.isLogInfos(sim.logger)
ModiaMath.Logging.isLogWarnings(sim::SimulationState) = ModiaMath.Logging.isLogWarnings(sim.logger)
ModiaMath.Logging.isLogEvents(sim::SimulationState) = ModiaMath.Logging.isLogEvents(sim.logger)



"""
    ModiaMath.setNominal!(sim::ModiaMath.SimulationState, x_nominal::Vector{Float64})

At initialization provide `x_nominal, the nominal values of the x-vector, and 
store them in `sim`, the simulation state. These values are used to compute the
absolute tolerances of the x-vector for the integrator.

# Example

```julia
function getModelResidues!(m::Model, t, x, derx, r, w)
   sim = m.simulationState
   if ModiaMath.isInitial(sim)
      ModiaMath.setNominal!(sim, [1.0, 1e-5, 1e8])   # if length(x)=3
   end
   ...
end
```
"""
function setNominal!(sim::SimulationState, x_nominal::Vector{Float64})
    @assert(length(x_nominal) == sim.nx)
    @assert(x_nominal[ indmin(x_nominal) ] > 0.0)
    for i in eachindex(x_nominal)
        sim.x_nominal[i] = x_nominal[i]
    end
end


"""
    ModiaMath.setNextEvent!(sim, nextEventTime; integratoToEvent=true, restart = ModiaMath.Restart)

At an event instant (isEvent(sim) == true) trigger the next time event. 

# Arguments

- `sim::ModiaMath.SimulationState`: (Internal) simulation state provided by simulation engine.
- `nextEventTime::Float64`: Time instant of the next time event.
- `integrateToEvent::Bool`: If true, the integrator integrates exactly to this event.
  If false, the integrator might integrate beyond `nextEventTime` (so the step size of the
  integrator is not influenced by this event). This option is only useful, if information is inquired
  at the event and restart=ModiaMath.NoRestart is used. The benefit is that the integrator is 
  practically not influenced by this event.
- `restart::ModiaMath.`[`EventRestart`](@ref): Restart action after the event at which
  `setNextEvent!` was called.
"""
setNextEvent!(sim::SimulationState, nextEventTime::Float64;
              integrateToEvent::Bool=true, restart::EventRestart=Restart) =
              setNextEvent!(sim.eventHandler, nextEventTime; integrateToEvent=integrateToEvent, restart=restart)



"""
    ModiaMath.positive!(sim, nr, crossing, crossingAsString; restart = ModiaMath.Restart)

Return `crossing > 0` such that a state event is triggered whenever `crossing > 0` changes its value.

Note, a small hysteresis is added to the crossing function during continuous integration, 
in order to reduce issues during event iterations due to small numerical errors.
However, at an event instant, `crossing > 0` is returned without any hysteresis.


# Arguments

- `sim::ModiaMath.SimulationState`: (Internal) simulation state provided by simulation engine.
- `nr::Int`: (> 0 required) Every call of `positive!(..)` must be identified by a unique value of `nr`
        This value is used as index in a vector that holds the internal memory for `crossing > 0`.
- `crossing::Float64`: Zero crossing function.
- `crossingAsString::String`: `crossing` as string representation. This string is used for log messages.
- `restart::ModiaMath.`[`EventRestart`](@ref): Restart action after the `crossing > 0` event occurred.


# Example

```julia
function computeVariables!(m::Model, sim::SimulationState)
   ...
   # f = s > 0 ? fMax : 0.0
   m.sPos.value = ModiaMath.positive!(sim, 1, m.s.value, "s")
   m.f.value    = m.sPos.value ? m.fMax.value : 0.0
   ...
end
```
"""
positive!(sim::SimulationState, nr::Int, crossing::Float64, crossingAsString::String; restart::EventRestart=Restart) =
              positive!(sim.eventHandler, nr, crossing, crossingAsString; restart=restart)



"""
    ModiaMath.negative!(sim, nr, crossing, crossingAsString; restart = ModiaMath.Restart)

Return `crossing < 0` such that a state event is triggered whenever `crossing < 0` changes its value.
"""
negative!(sim::SimulationState, nr::Int, crossing::Float64, crossingAsString::String; restart::EventRestart=Restart) =
              negative!(sim.eventHandler, nr, crossing, crossingAsString; restart=restart)      

  

"""
    ModiaMath.change!(sim, nr, crossing, crossingAsString; restart = ModiaMath.Restart)

Trigger an event, whenever `crossing > 0` changes from `false` to `true` or from `true` to `false`.
The function returns always `false`.
"""      
change!(sim::SimulationState, nr::Int, crossing::Float64, crossingAsString::String; restart::EventRestart=Restart) =
              change!(getDaeInfo(model).eventHandler, nr, crossing, crossingAsString; restart=restart)


"""
    ModiaMath.edge!(sim, nr, crossing, crossingAsString; restart = ModiaMath.Restart)

Trigger an event, whenever `crossing > 0` switches from `false` to `true`.
The function returns always `false`.
"""
edge!(sim::SimulationState, nr::Int, crossing::Float64, crossingAsString::String; restart::EventRestart=Restart) =
              edge!(sim.eventHandler, nr, crossing, crossingAsString; restart=restart)



ModiaMath.Logging.logOn!(sim::SimulationState)                     = ModiaMath.Logging.logOn!(sim.logger)  
ModiaMath.Logging.logOn!(model::ModiaMath.AbstractSimulationModel) = ModiaMath.Logging.logOn!(model.simulationState.logger) 

ModiaMath.Logging.logOff!(sim::SimulationState)                     = ModiaMath.Logging.logOff!(sim.logger)
ModiaMath.Logging.logOff!(model::ModiaMath.AbstractSimulationModel) = ModiaMath.Logging.logOff!(model.simulationState.logger)


"""
    ModiaMath.setLogCategories!(m::ModiaMath.[AbstractSimulationModel|SimulationState], categories; reinit=true)

Set log categories as vector of symbols, e.g. setLogCategories!(logger, [:LogProgess]).
Supported categories:
- `:LogStatistics`, print statistics information at end of simulation
- `:LogProgress`, print progress information during simulation
- `:LogInfos`, print information messages of the model
- `:LogWarnings`, print warning messages of the model
- `:LogEvents`, log events of the model

If option reinit=true, all previous set categories are reinitialized to be no longer present.
If reinit=false, previously set categories are not changed.
"""
ModiaMath.Logging.setLogCategories!(sim::SimulationState, categories::Vector{Symbol}; reinit::Bool=true) = 
                  ModiaMath.Logging.setLogCategories!(sim.logger, categories; reinit=reinit)
ModiaMath.Logging.setLogCategories!(model::ModiaMath.AbstractSimulationModel, categories::Vector{Symbol}; reinit::Bool=true) = 
                  ModiaMath.Logging.setLogCategories!(model.simulationState.logger, categories; reinit=reinit)


# Functions that can be called after a model is instantiated
getSimulationResult(model::ModiaMath.AbstractSimulationModel) = model.simulationState.result
