# License for this file: MIT (expat)
# Copyright 2020, DLR Institute of System Dynamics and Control

"""
    module ModiaTest_ode_pendulum

ODE-model of a mass point attached via a rod to a revolute joint with 2 states.

This test is used to evalute how to use DifferentialEquations.jl in ModiaMath
utilizing the new Modia structure that is currently under development.
"""
module ModiaTest_ode_pendulum

import DifferentialEquations
using PyPlot        # For plotting

using StaticArrays  # Test of StaticArrays with parameterized data
using DoubleFloats  # Test with higher precision (128 bit instead of 64 bit floating point)
using Measurements  # Test with uncertain variables (nominal value + uncertainty)


#= Number types supported by DifferentialEquations.jl: https://tutorials.juliadiffeq.org/html/type_handling/01-number_types.html

using DoubleFloat
DoubleFloat{Float64}: 128 bit (34 digits)

using Measurements
Measurement{FLoat64}: 64bit ± 64bit

using ArbNumerics
setworkingprecision(ArbFloat, bits=500)
ArbFloat{500}: 500 bit
ArbReal{500} : 500 bit with precision interval
               Can be set explicitly
               setball(midpoint::ArbFloat, radius::ArbFloat)
               setinterval(lobound::ArbFloat, hibound::ArbFloat)
               (midpoint, radius) = ball(x::ArbReal)
               (lobound, hibound) = interval(x::ArbReal)

using BigFloat   // ArbNumerics is faster, but BigFLoat is supported by more Julia functions
=#


#---------------------------------------- ModiaMath definitions

"""
    AbstractSimulationModel{StateType,RealType}

Supertype of all simulation models

- StateType: Type of the state vector and of its derivative:
                x::Vector{StateType}
                der_x::Vector{StateType}
             and all variables that are used to compute der_x from x.
             Parameters can be either of StateType (e.g. to support parameters with uncertainty)
             or RealType (if not to be included in sensitivity analysis)

- RealType : Type of floating point variables that are not of type StateType, such as:
             time, tolerance, startTime, StopTime, interval, nominal, ...
"""
abstract type AbstractSimulationModel{StateType,RealType} end


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
    Result                   # constructor for result data structure
    addToResult!::Function

    # Problem variables
    tolerance::RealType           # Relative integration tolerance
    startTime::RealType           # Start time of the simulation
    stopTime::RealType            # Stop time of the simulation
    interval::RealType            # Interval of the simulation

    # Internal state
    storeResult::Bool

    # Reference to simulation model
    simulationModel

    # Simulation result
    result

    SimulationState{StateType,RealType}(modelName, x_start,
                                        setStates!, getDerivatives!, setStatesAndDerivatives!, getResiduals!,
                                        evaluate!,  outputs!, Result, addToResult!;
                                        tolerance=1e-4, startTime=0.0, stopTime=1.0, interval=0.02, dae=true) where {StateType, RealType} =
         new(modelName, x_start, similar(x_start), similar(x_start), dae, setStates!, getDerivatives!, setStatesAndDerivatives!, getResiduals!,
             evaluate!, outputs!, Result, addToResult!, tolerance, startTime, stopTime, interval, false)
end


#
# x, der_x, residuals need not to be Vector{..}, e.g. Sundials.IDA provides N_Vector
#

"derivatives!: Called by ODE integrator for ODE model"
function derivatives!(der_x::AbstractVector, x::AbstractVector, simulationModel, t)::Nothing
    sim::SimulationState = simulationModel.simulationState
    sim.setStates!(simulationModel, x)
    sim.evaluate!(simulationModel)
    sim.getDerivatives!(simulationModel, der_x)
    return nothing
end


"ODEresiduals!: Called by DAE integrator for ODE model"
function ODEresiduals!(residuals::AbstractVector, der_x::AbstractVector, x::AbstractVector, simulationModel, t)::Nothing
    sim::SimulationState = simulationModel.simulationState
    derivatives!(sim.der_x, x, simulationModel, t)
    for i = 1:length(der_x)
        residuals[i] = sim.der_x[i] - der_x[i]
    end
    return nothing
end


"ODEoutputs!: Called by ODE or DAE integrator for ODE model at every communication point"
function ODEoutputs!(x::AbstractVector, t, integrator)::Nothing
    simulationModel = integrator.p
    sim::SimulationState = simulationModel.simulationState
    sim.setStates!(simulationModel, x)
    sim.storeResult = true
    sim.evaluate!(simulationModel)
    sim.storeResult = false
    sim.addToResult!(sim.result, t, simulationModel)
    return nothing
end


function setStatesAndDerivativesDummy!(simulationModel, x::AbstractVector, der_x::AbstractVector)::Nothing
    error("setStatesAndDerivativesDummy! should never be called")
    return nothing
end


"SimulationState: Constructor for ODE Simulation State"
SimulationState{StateType,RealType}(modelName, x_start, setStates!, getDerivatives!, evaluate!, Result, addToResult!;
                                    tolerance=1e-4, startTime=0.0, stopTime=1.0, interval=0.02) where {StateType, RealType} =
    SimulationState{StateType,RealType}(modelName, x_start, setStates!, getDerivatives!, setStatesAndDerivativesDummy!, ODEresiduals!,
                                        evaluate!,  ODEoutputs!, Result, addToResult!; dae=false,
                                        tolerance=tolerance, startTime=startTime, stopTime=stopTime, interval=interval)


"solve!: Return DifferentialEquations.solve(...)"
function solve!(simulationModel::AbstractSimulationModel{StateType,RealType}, integrator; tolerance=NaN, startTime=NaN, stopTime=NaN, interval=NaN) where {StateType, RealType}
    sim::SimulationState{StateType,RealType} = simulationModel.simulationState

    # Construct a new result instance
    sim.result = sim.Result{StateType,RealType}(simulationModel)

    # Change tolerance, startTime, stopTime, interval if required
    if !isnan(tolerance); @assert(tolerance > 0); sim.tolerance = RealType(tolerance); end
    if !isnan(startTime); sim.startTime = RealType(startTime); end
    if !isnan(stopTime) ; sim.stopTime  = RealType(stopTime) ; end
    if !isnan(interval) ; @assert(interval > 0); sim.interval = RealType(interval) ; end

    # Define problem based on integrator and model type
    if typeof(integrator) <: DifferentialEquations.IDA
        # DAE integrator
        if !sim.dae
            # ODE model
            derivatives!(sim.der_x_start, sim.x_start, simulationModel, sim.startTime)
            problem = DifferentialEquations.DAEProblem(ODEresiduals!, sim.der_x_start, sim.x_start, (sim.startTime, sim.stopTime), simulationModel)
        else
            # DAE model
            error("DAE model not yet implemented")
        end
    else
        # Hope that it is an ODE integrator
        problem = DifferentialEquations.ODEProblem(derivatives!, sim.x_start, (sim.startTime, sim.stopTime), simulationModel)
    end

    # Callbacks
    tspan = sim.startTime:sim.interval:sim.stopTime
    cb1 = DifferentialEquations.FunctionCallingCallback(sim.outputs!, funcat=tspan)

    # Simulation
    println("... Simulation started for model ", sim.modelName, " using integrator ", typeof(integrator), " with StateType = ", string(StateType) )
    solution = DifferentialEquations.solve(problem, integrator, reltol=sim.tolerance, saveat=sim.stopTime-sim.startTime, callback=cb1)
end



#------------------------------------------------------------------------------------------------
#---------------------------------------- Provided by Modia (rough sketch as needed for testing)
#------------------------------------------------------------------------------------------------

#--------------------------------------- Generic for all models

mutable struct Variable
    typ
    start
    fixed::Union{Bool,Nothing}
    info::AbstractString

    Variable(;typ=Real, start=nothing, fixed=nothing, info="") = new(typ, start, fixed, info)
end
const VariableDict = Dict{Symbol, Variable}


mutable struct ModelResult{StateType,RealType}
    t::Vector{RealType}
    model::AbstractVector
    variables::VariableDict

    ModelResult{StateType,RealType}(m::AbstractSimulationModel) where {StateType, RealType} =
          new(RealType[], typeof(m.model)[], m.variables)
end


function addToResult!(result::ModelResult{StateType,RealType}, t, m::AbstractSimulationModel{StateType,RealType})::Nothing where {StateType, RealType}
    push!(result.t, deepcopy(t))
    push!(result.model, deepcopy(m.model))
    return nothing;
end



#---------------------------------------- Specific to a particular model (here: Pendulum)

mutable struct Model{StateType,RealType}
    # Parameters
    L::StateType
    m::StateType
    d::StateType
    g::StateType

    # Variables during integration
    phi::StateType
    w::StateType
    der_w::StateType

    # Variables at communication points
    r::SVector{2,StateType}   # absolute position of end point

    function Model{StateType,RealType}(;L=1.0, m=1.0, d=0.1, g=9.81) where {StateType, RealType}
        @assert(L > 0)
        @assert(m > 0)
        @assert(d >= 0)
        new(L, m, d, g, 0, 0, 0, SVector{2,StateType}(0,0))
    end
end


mutable struct PendulumGeneric{StateType,RealType} <: AbstractSimulationModel{StateType,RealType}
    # Required:
    simulationState::SimulationState{StateType, RealType}

    # Defined by Modia
    model::Model{StateType,RealType}   # Data structure of model variable values
    variables::VariableDict            # Data structure of model variable definition
end
const Pendulum = PendulumGeneric{Float64, Float64}


# The following definition tells DifferentialEquations.jl that a SimulationModel is mutable (otherwise stack overflow error for two calls of solve!)
using ArrayInterface
ArrayInterface.ismutable(::Type{PendulumGeneric}) = true



function ModelVariables()::VariableDict
    var = VariableDict()

    # Parameters
    var[:L] = Variable(info="Angle")
    var[:m] = Variable(info="Mass")
    var[:d] = Variable(info="Damping factor")
    var[:g] = Variable(info="Gravity acceleration")

    # Variables during integration
    var[:phi]   = Variable(start=π/2, info="Angle")
    var[:w]     = Variable(start=0.0, info="Angular Velocity (= der(phi))")
    var[:der_w] = Variable(info="Angular Acceleration (= der(w))")

    # Variables at communication points
    var[:r]     = Variable(typ=SVector{2,Real}, info="Absolute position of pendulum tip")

    return var
end


function InitialStates(model::Model{StateType,RealType}, var::VariableDict)::Vector{StateType} where {StateType, RealType}
    x_start = [StateType(var[:phi].start),
               StateType(var[:w].start)]
    return x_start
end

function setStates!(obj::PendulumGeneric{StateType,RealType}, x::AbstractVector)::Nothing where {StateType, RealType}
    m = obj.model
    m.phi = x[1]
    m.w   = x[2]
    return nothing
end

function getDerivatives!(obj::PendulumGeneric{StateType,RealType}, der_x::AbstractVector)::Nothing where {StateType, RealType}
    m = obj.model
    der_x[1] = m.w
    der_x[2] = m.der_w
    return nothing
end

function evaluate!(obj::PendulumGeneric{StateType,RealType})::Nothing where {StateType, RealType}
    m   = obj.model
    sim = obj.simulationState
    m.der_w = (-m.m * m.g * m.L * sin(m.phi) - m.d * m.w) / (m.m * m.L^2)

    if sim.storeResult
        m.r = SVector{2,StateType}(m.L*sin(m.phi), -m.L*cos(m.phi))
    end
    return nothing
end

"Pendulum: Constructor for Pendulum SimulationModel"
function PendulumGeneric{StateType,RealType}(;kwargs...)  where {StateType, RealType}
    model     = Model{StateType,RealType}(;kwargs...)
    variables = ModelVariables()
    x_start   = InitialStates(model, variables)

    simulationState = SimulationState{StateType,RealType}("Pendulum", x_start, setStates!, getDerivatives!, evaluate!, ModelResult, addToResult!)

    simulationModel = PendulumGeneric{StateType,RealType}(simulationState, model, variables)
    simulationState.simulationModel = simulationModel
    return simulationModel
end




#------------------------------------------------------------------------------------------------
#---------------------------------------- Simulation of the pendulum model
#------------------------------------------------------------------------------------------------

function standardPlot(solution, fig1, fig2)
    result = solution.prob.p.simulationState.result
    StateType = string( typeof( solution.prob.p.simulationState.x_start[1] ) )
    figure(fig1)
    clf()
    t   = result.t
    phi = getfield.(result.model, :phi)
    plot(t, phi, label="\$\\varphi\$")
    grid(true)
    legend()
    title("SimplePendulum with " * StateType)

    figure(fig2)
    clf()
    r  = getfield.(result.model, :r)
    r1 = getindex.(r,1)
    r2 = getindex.(r,2)
    plot(t, r1, label="r[1]")
    plot(t, r2, label="r[2]")
    grid(true)
    legend()
    title("SimplePendulum with " * StateType)
end


# ODE integrator with Float64
pendulum = Pendulum(m=1.1)
solution = solve!(pendulum, DifferentialEquations.Tsit5(), stopTime=7.0, tolerance=1e-4)
standardPlot(solution, 1, 2)

# DAE integrator with FLoat64
solution = solve!(pendulum, DifferentialEquations.IDA(), stopTime=7.0, tolerance=1e-4)
standardPlot(solution, 3, 4)

# ODE integrator with DoubleFloats
pendulum = PendulumGeneric{DoubleFloat{Float64},Float64}(m=1.1)
solution = solve!(pendulum, DifferentialEquations.Tsit5(), stopTime=7.0, tolerance=1e-4)
standardPlot(solution, 5, 6)

# ODE integrator with Measurements
pendulum = PendulumGeneric{Measurement{Float64},Float64}(L=1.0±0.1, m=1.0±0.1)
solution = solve!(pendulum, DifferentialEquations.Tsit5(), stopTime=7.0, tolerance=1e-4)
result   = solution.prob.p.simulationState.result
StateType = string( typeof( solution.prob.p.simulationState.x_start[1] ) )

t   = result.t
phi = getfield.(result.model, :phi)
phi_mean = Measurements.value.(phi)
phi_u    = Measurements.uncertainty.(phi)
phi_min  = phi_mean + phi_u
phi_max  = phi_mean - phi_u

figure(7)
clf()
plot(t, phi_max , label="\$\\varphi_{max}\$")
plot(t, phi_min , label="\$\\varphi_{min}\$" )
plot(t, phi_mean, label="\$\\varphi_{mean}\$" )
grid(true)
legend()
title("SimplePendulum with Measurement{Float64}")


end