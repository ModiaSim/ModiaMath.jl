# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module ModiaMath.DAE

Interface between the [`ModiaMath.SimulationEngine`](@ref) and the index 1 DAE model.

The following functions can be called in the DAE model to inquire
**information about the simulation state**:

- [`ModiaMath.getTime`](@ref)
- [`ModiaMath.getStartTime`](@ref)
- [`ModiaMath.getStopTime`](@ref)
- [`ModiaMath.getTolerance`](@ref)
- [`ModiaMath.isInitial`](@ref)
- [`ModiaMath.isTerminal`](@ref)
- [`ModiaMath.isEvent`](@ref)
- [`ModiaMath.isZeroCrossing`](@ref)
- [`ModiaMath.isAfterSimulationStart`](@ref)
- [`ModiaMath.isStoreResult`](@ref)
- [`ModiaMath.isLogInfos`](@ref)
- [`ModiaMath.isLogWarnings`](@ref)
- [`ModiaMath.isLogEvents`](@ref)

The following functions can be called in the DAE model to set
**properties in the simulation engine**:

- [`ModiaMath.setNominal!`](@ref)
- [`ModiaMath.setNextEvent!`](@ref)
- [`ModiaMath.positive!`](@ref)
- [`ModiaMath.negative!`](@ref)
- [`ModiaMath.change!`](@ref)
- [`ModiaMath.edge!`](@ref)

The following functions can be either called in the DAE model
or they can be called on a simulation model
(before or after `ModiaMath.simulate!(simulationModel, ...)` is called).

- [`ModiaMath.logOn!`](@ref)
- [`ModiaMath.logOff!`](@ref)
- [`ModiaMath.setLogCategories!`](@ref)


# Main developer

[Martin Otter](https://rmc.dlr.de/sr/de/staff/martin.otter/),
[DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)
"""
module DAE

export EventRestart, NoRestart, Restart, FullRestart, Terminate
export getVariableName, getResultNames

export getTime, getStartTime, getStopTime, getTolerance
export isInitial, isTerminal, isEvent, isZeroCrossing, isAfterSimulationStart, isStoreResult, setNominal!
export setNextEvent!, positive!, negative!, change!, edge!
export getSimulationResult

export SimulationState

# export InitInfo, EventInfo, reset!
# export EventHandler
# export getEventResidues!
# export reinitialize!, eventIteration!, initialize!, getResidues!
# export computeAndStoreResult!, terminate!, processEvent!, getEventIndicators

import ModiaMath
using LinearAlgebra

"""
    @enum EventRestart NoRestart Restart FullRestart Terminate

Define how to continue or restart integration after an event. Usually, `Restart`
should be used. Only in special cases, the other flags are useful.

- `NoRestart`, continue integration without restarting the integrator
- `Restart`, restart integrator
- `FullRestart`, restart integrator and optionally exchange simulationState (so dimensions may change)
- `Terminate`, terminate integration
"""
@enum EventRestart NoRestart Restart FullRestart Terminate

@enum VariableCategory Category_X Category_DERX Category_W


# include code
include("events.jl")
include("returnStructs.jl")
include("simulationState.jl")
include("functionsForUserModels.jl")

end
