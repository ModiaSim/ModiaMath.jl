# License for this file: MIT (expat)
# Copyright 2020, DLR Institute of System Dynamics and Control

"""
    module ModiaMathTest

Evaluate the basic structure of the ModiaMath simulator supporting DifferentialEquations.jl
"""
module ModiaMathTest

import DifferentialEquations
export SimulationState, solve!


#---------------------------------------- ModiaMath definitions

#=
Used Types

- StateType: Type of the state vector, of its derivative and of time:
                x::Vector{StateType}
                der_x::Vector{StateType}
				time::StateType
             and all variables that are used to compute der_x from x and time.
             Parameters can be either of StateType (e.g. to support parameters with uncertainty)
             or RealType (if not to be included in sensitivity analysis)

- RealType : Type of floating point variables that are not of type StateType, such as:
             tolerance, startTime, StopTime, interval, nominal, ...
=#




mutable struct SimulationState{StateType,RealType}
    modelName::AbstractString

    x_start::Vector{StateType}
    der_x_start::Vector{StateType}

    # Cache
    der_x::Vector{StateType}

    # Type of model
    dae::Bool    # = true: DAE; =false: ODE

    # Functions for ODE models
    setStates!::Function
    getDerivatives!::Function

    # Functions for DAE models
    setStatesAndDerivatives!::Function

    # Functions for ODE and DAE models
    getResiduals!::Function
    evaluate!::Function
    outputs!::Function
    Result                        # constructor for result data structure
    addToResult!::Function

    # Problem variables
    tolerance::RealType           # Relative integration tolerance
    startTime::RealType           # Start time of the simulation
    stopTime::RealType            # Stop time of the simulation
    interval::RealType            # Interval of the simulation

    # Time varying variables
    time::StateType
    storeResult::Bool

    # Simulation result
    result::Any                   # Instance of result

    SimulationState{StateType,RealType}(modelName, x_start,
                                        setStates!, getDerivatives!, setStatesAndDerivatives!, getResiduals!,
                                        evaluate!,  outputs!, Result, addToResult!;
                                        tolerance=1e-4, startTime=0.0, stopTime=1.0, interval=NaN, dae=true) where {StateType, RealType} =
         new(modelName, x_start, similar(x_start), similar(x_start), dae, setStates!, getDerivatives!, setStatesAndDerivatives!, getResiduals!,
             evaluate!, outputs!, Result, addToResult!, tolerance, startTime, stopTime, interval, RealType(0), false)
end


#
# x, der_x, residuals need not to be Vector{..}, e.g. Sundials.IDA provides N_Vector
#

"derivatives!: Called by ODE integrator for ODE model"
function derivatives!(der_x, x, simulationModel, t)::Nothing
    sim::SimulationState = simulationModel.simulationState
    sim.time = t
    sim.setStates!(simulationModel, x)
    sim.evaluate!(simulationModel)
    sim.getDerivatives!(simulationModel, der_x)
    return nothing
end


"ODEresiduals!: Called by DAE integrator for ODE model"
function ODEresiduals!(residuals, der_x, x, simulationModel, t)::Nothing
    sim::SimulationState = simulationModel.simulationState
    sim.time = t
    derivatives!(sim.der_x, x, simulationModel, t)
    for i = 1:length(der_x)
        residuals[i] = sim.der_x[i] - der_x[i]
    end
    return nothing
end


"ODEoutputs!: Called by ODE or DAE integrator for ODE model at every communication point"
function ODEoutputs!(x, t, integrator)::Nothing
    simulationModel = integrator.p
    sim::SimulationState = simulationModel.simulationState
    sim.time = t
    sim.setStates!(simulationModel, x)
    sim.storeResult = true
    sim.evaluate!(simulationModel)
    sim.storeResult = false
    sim.addToResult!(sim.result, t, simulationModel)
    return nothing
end

"initialize!: Called for DAE model to initialize"
function initialize!(simulationModel, sim::SimulationState{StateType,RealType})::Nothing where {StateType,RealType}
    if sim.dae
		error("... DAE initialization not yet implemented")
	else
		println("    ODE initialization function called")
		derivatives!(sim.der_x_start, sim.x_start, simulationModel, StateType(sim.startTime))
	end
    return nothing
end

function setStatesAndDerivativesDummy!(simulationModel, x, der_x)::Nothing
    error("setStatesAndDerivativesDummy! should never be called")
    return nothing
end


"SimulationState: Constructor for ODE Simulation State"
SimulationState{StateType,RealType}(modelName, x_start, setStates!, getDerivatives!, evaluate!, Result, addToResult!;
                                    tolerance=1e-4, startTime=0.0, stopTime=1.0, interval=NaN) where {StateType, RealType} =
    SimulationState{StateType,RealType}(modelName, x_start, setStates!, getDerivatives!, setStatesAndDerivativesDummy!, ODEresiduals!,
                                        evaluate!,  ODEoutputs!, Result, addToResult!; dae=false,
                                        tolerance=tolerance, startTime=startTime, stopTime=stopTime, interval=interval)

"solve!: Return DifferentialEquations.solve(...)"
function solve!(simulationModel, integrator; tolerance=NaN, startTime=NaN, stopTime=NaN, interval=NaN)
    sim::SimulationState = simulationModel.simulationState
    RealType  = typeof(sim.tolerance)
    StateType = typeof(sim.time)
    println("... Simulation started for model ", sim.modelName, " using integrator ", typeof(integrator), " with StateType = ", string(StateType) )

    # Construct a new result instance
    sim.result = sim.Result(simulationModel)

    # Change tolerance, startTime, stopTime, interval if required
    if !isnan(tolerance); @assert(tolerance > 0); sim.tolerance = RealType(tolerance); end
    if !isnan(startTime); sim.startTime = RealType(startTime); end
    if !isnan(stopTime) ; sim.stopTime  = RealType(stopTime) ; end
    if isnan(interval)
		if isnan(sim.interval)
			sim.interval = abs(sim.stopTime - sim.startTime)/500
		end
	else
		@assert(interval > 0)
		sim.interval = RealType(interval)
    end

    # Initialize model
    sim.time = StateType(sim.startTime)
    initialize!(simulationModel, sim)

    # Define problem based on integrator and model type
    if typeof(integrator) <: DifferentialEquations.IDA
        # DAE integrator
        if !sim.dae
            # ODE model
            problem = DifferentialEquations.DAEProblem(ODEresiduals!, sim.der_x_start, sim.x_start, (StateType(sim.startTime), sim.stopTime), simulationModel)
        else
            # DAE model
            error("DAE model not yet implemented")
        end
    else
        # Hope that it is an ODE integrator
        problem = DifferentialEquations.ODEProblem(derivatives!, sim.x_start, (StateType(sim.startTime), sim.stopTime), simulationModel)
    end

    # Callbacks
    tspan = sim.startTime:sim.interval:sim.stopTime
    cb1 = DifferentialEquations.FunctionCallingCallback(sim.outputs!, funcat=tspan)

    # Simulation
    solution = DifferentialEquations.solve(problem, integrator, reltol=sim.tolerance, saveat=sim.stopTime-sim.startTime, callback=cb1)
end



end