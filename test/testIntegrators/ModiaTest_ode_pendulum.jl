# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module ModiaTest_ode_pendulum

ODE-model of a mass point attached via a rod to a revolute joint with 2 states.
"""
module ModiaTest_ode_pendulum

import DifferentialEquations
using  StaticArrays

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
    AbstractModel{StateType<:Real,RealType<:Real}

Supertype of all models
"""
abstract type AbstractModel{StateType<:Real,RealType<:Real} end


mutable struct SimulationState{StateType<:Real, RealType<:Real}
    x_start::Vector{StateType}
    der_x_start::Vector{StateType}

    # Cache
    der_x::Vector{StateType}

    # Functions for ODE models
    setStates!::Function
    getDerivatives!::Function

    # Functions for DAE models
    setStatesAndDerivatives!::Function

    # Functions for ODE and DAE models
    getResiduals!::Function
    evaluate!::Function
    outputs!::Function
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
    ode::Bool    # = true: ODE; =false: DAE

    SimulationState{StateType,RealType}(x_start,
                                        setStates!, getDerivatives!, setStatesAndDerivatives!, getResiduals!,
                                        evaluate!,  outputs!, addToResult!,
                                        tolerance, startTime, stopTime, interval, ode) where {StateType<:Real, RealType<:Real} =
         new(x_start, similar(x_start), similar(x_start), setStates!, getDerivatives!, setStatesAndDerivatives!, getResiduals!,
             evaluate!, outputs!, addToResult!, tolerance, startTime, stopTime, interval, false, nothing, ode)
end


mutable struct SimulationModel{StateType,RealType}
    modelName::AbstractString                  # Name of top level model instance
    model::AbstractModel{StateType,RealType}   # Data structure of model variable values
    variables                                  # Data structure of model variable definition
    result                                     # Data structure of result
    simulationState::SimulationState{StateType, RealType}

    function SimulationModel{StateType,RealType}(modelName, model::AbstractModel{StateType,RealType}, variables, result,
                             simulationState::SimulationState{StateType,RealType}) where {StateType<:Real, RealType<:Real}
        obj = new{StateType,RealType}(modelName, model, variables, result, simulationState)
        simulationState.simulationModel = obj
        return obj
    end
end

# The following definition tells DifferentialEquations.jl that a SimulationModel is mutable (otherwise stack overflow error)
using ArrayInterface
ArrayInterface.ismutable(::Type{SimulationModel}) = true


"derivatives!: Called by ODE integrator for ODE model"
function derivatives!(der_x, x, simulationModel::SimulationModel, t)::Nothing
    sim = simulationModel.simulationState
    sim.setStates!(simulationModel, x)
    sim.evaluate!(simulationModel)
    sim.getDerivatives!(simulationModel, der_x)
    return nothing
end


"ODEresiduals!: Called by DAE integrator for ODE model"
function ODEresiduals!(residuals, der_x, x, simulationModel::SimulationModel, t)::Nothing
    sim = simulationModel.simulationState
    derivatives!(sim.der_x, x, simulationModel, t)
    for i = 1:length(der_x)
        residuals[i] = sim.der_x[i] - der_x[i]
    end
    return nothing
end


"ODEoutputs!: Called by ODE or DAE integrator for ODE model at every communication point"
function ODEoutputs!(x, t, integrator)::Nothing
    simulationModel = integrator.p
    sim = simulationModel.simulationState
    sim.setStates!(simulationModel, x)
    sim.storeResult = true
    sim.evaluate!(simulationModel)
    sim.storeResult = false
    sim.addToResult!(simulationModel.result, t, simulationModel.model)
    return nothing
end


function setStatesAndDerivativesDummy!(simModel::SimulationModel{StateType,RealType}, x::Vector{StateType}, der_x::Vector{StateType})::Nothing where {StateType<:Real, RealType<:Real}
    error("setStatesAndDerivativesDummy! should never be called")
    return nothing
end


"ODESimulationModel: Constructor for ODE SimulationModel"
function ODESimulationModel(modelName, model::AbstractModel{StateType,RealType}, variables, result, x_start::Vector{StateType}, setStates!, getDerivatives!, evaluate!, addToResult!;
                            tolerance=1/10^4, startTime=0.0, stopTime=1.0, interval=NaN) where {StateType<:Real, RealType<:Real}

    @assert(tolerance > 0)
    @assert(stopTime >= startTime)
    tolerance2 = RealType(tolerance)
    startTime2 = RealType(startTime)
    stopTime2  = RealType(stopTime)
    interval2  = isnan(interval) ? (stopTime2 - startTime2)/RealType(500) : RealType(interval)

    simulationState = SimulationState{StateType,RealType}(x_start, setStates!, getDerivatives!, setStatesAndDerivativesDummy!,
                                      ODEresiduals!, evaluate!, ODEoutputs!, addToResult!,
                                      tolerance2, startTime2, stopTime2, interval2, true)

    simulationModel = SimulationModel{StateType,RealType}(modelName, model, variables, result, simulationState)
end


"solve!: Return DifferentialEquations.solve(...)"
function solve!(simulationModel::SimulationModel{StateType,RealType}, integrator; tolerance=NaN, startTime=NaN, stopTime=NaN, interval=NaN) where {StateType<:Real, RealType<:Real}
    sim = simulationModel.simulationState

    # Change tolerance, startTime, stopTime, interval if required
    if !isnan(tolerance); @assert(tolerance > 0); sim.tolerance = RealType(tolerance); end
    if !isnan(interval) ; @assert(interval  > 0); sim.interval  = RealType(interval) ; end
    if !isnan(startTime); sim.startTime = RealType(startTime); end
    if !isnan(stopTime) ; sim.stopTime  = RealType(stopTime) ; @assert(sim.stopTime >= sim.StartTime); end

    # Define problem based on integrator and model type
    if typeof(integrator) <: DifferentialEquations.IDA
        # DAE integrator
        if sim.ode
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
    println("... Simulation started for model ", simulationModel.modelName, " using integrator ", typeof(integrator), " with StateType = ", string(StateType) )
    solution = DifferentialEquations.solve(problem, integrator, reltol=sim.tolerance, saveat=sim.stopTime-sim.startTime, callback=cb1)
end



#---------------------------------------- General definitions for models (e.g. provided by Modia)
mutable struct Variable
    typ
    start
    fixed::Union{Bool,Nothing}
    info::AbstractString

    Variable(;typ=Real, start=nothing, fixed=nothing, info="") = new(typ, start, fixed, info)
end

const VariableDict = Dict{Symbol, Variable}



#---------------------------------------- Specific to a particular model (e.g. provided by Modia)

mutable struct Model{StateType,RealType} <: AbstractModel{StateType,RealType}
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

    function Model{StateType,RealType}(;L=1.0, m=1.0, d=0.1, g=9.81) where {StateType<:Real, RealType<:Real}
        @assert(L > 0)
        @assert(m > 0)
        @assert(d >= 0)
        new(L, m, d, g, 0, 0, 0, SVector{2,StateType}(0,0))
    end
end


#=
mutable struct ModelResult{StateType<:Real, RealType<:Real}
    t::Vector{RealType}
    model::Vector{Model{StateType,RealType}}
    ModelResult{StateType,RealType}() where {StateType<:Real, RealType<:Real} = new(RealType[], Model{StateType,RealType}[])
end
=#

mutable struct ModelResult{StateType<:Real, RealType<:Real}
    t::Vector{RealType}
    model::Vector{Model{StateType,RealType}}
    variables::VariableDict

    ModelResult(m::Model{StateType,RealType}, var::VariableDict) where {StateType<:Real, RealType<:Real} = new{StateType,RealType}(RealType[], Model{StateType,RealType}[], var)
end


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



function InitialStates(m::Model{StateType,RealType}, var::VariableDict)::Vector{StateType} where {StateType<:Real, RealType<:Real}
    x_start = [StateType(var[:phi].start),
               StateType(var[:w].start)]
    return x_start
end

function setStates!(simModel::SimulationModel{StateType,RealType}, x)::Nothing where {StateType<:Real, RealType<:Real}
    m = simModel.model
    m.phi = x[1]
    m.w   = x[2]
    return nothing
end


function getDerivatives!(simModel::SimulationModel{StateType,RealType}, der_x::Vector{StateType})::Nothing where {StateType<:Real, RealType<:Real}
    m = simModel.model
    der_x[1] = m.w
    der_x[2] = m.der_w
    return nothing
end

function evaluate!(simModel::SimulationModel{StateType,RealType})::Nothing where {StateType<:Real, RealType<:Real}
    m   = simModel.model
    sim = simModel.simulationState
    m.der_w = (-m.m * m.g * m.L * sin(m.phi) - m.d * m.w) / (m.m * m.L^2)

    if sim.storeResult
        m.r = SVector{2,StateType}(m.L*sin(m.phi), -m.L*cos(m.phi))
    end
    return nothing
end

function addToResult!(result::ModelResult{StateType,RealType}, t::RealType, model::Model{StateType,RealType})::Nothing where {StateType<:Real, RealType<:Real}
    push!(result.t, t)
    push!(result.model, deepcopy(model))
    return nothing;
end

using PyPlot

function standardPlot(solution, fig1, fig2)
    result = solution.prob.p.result
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

#------------------------------------ Simulation of a particular Model

# ODE integrator with Float64
model     = Model{Float64,Float64}(m=1.1)
variables = ModelVariables()
result    = ModelResult(model, variables)
x_start   = InitialStates(model, variables)
simulationModel = ODESimulationModel("PendulumODE", model, variables, result, x_start, setStates!, getDerivatives!, evaluate!, addToResult!, stopTime=7.0, tolerance=1e-6)
solution = solve!(simulationModel, DifferentialEquations.Tsit5())
standardPlot(solution, 1, 2)

# DAE integrator with FLoat64
solution  = solve!(simulationModel, DifferentialEquations.IDA())
standardPlot(solution, 3, 4)


# ODE integrator with DoubleFloats
using DoubleFloats
model     = Model{DoubleFloat{Float64},Float64}(m=1.1)
variables = ModelVariables()
result    = ModelResult(model, variables)
x_start   = InitialStates(model, variables)
simulationModel = ODESimulationModel("PendulumODE", model, variables, result, x_start, setStates!, getDerivatives!, evaluate!, addToResult!, stopTime=7.0, tolerance=1e-6)
solution = solve!(simulationModel, DifferentialEquations.Tsit5())
standardPlot(solution, 5, 6)


# ODE integrator with Measurements
using Measurements
model     = Model{Measurement{Float64},Float64}(L=1.0±0.1, m=1.0±0.1)
variables = ModelVariables()
result    = ModelResult(model, variables)
x_start   = InitialStates(model, variables)
simulationModel = ODESimulationModel("PendulumODE", model, variables, result, x_start, setStates!, getDerivatives!, evaluate!, addToResult!, stopTime=7.0, tolerance=1e-6)
solution = solve!(simulationModel, DifferentialEquations.Tsit5())


result = solution.prob.p.result
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