# License for this file: MIT (expat)
# Copyright 2020, DLR Institute of System Dynamics and Control

"""
    module ModiaMathTest_ode_pendulum

ODE-model of a mass point attached via a rod to a revolute joint with 2 states.

This test is used to evalute how to use DifferentialEquations.jl in ModiaMath
using different types: Float64,  DoubleFloat{Float64}, Measurement{Float64}
"""
module ModiaMathTest_ode_pendulum

import DifferentialEquations
using PyPlot        # For plotting

using StaticArrays  # Test of StaticArrays with parameterized data
using DoubleFloats  # Test with higher precision (128 bit instead of 64 bit floating point)
using Measurements  # Test with uncertain variables (nominal value + uncertainty)


include("ModiaMathTest.jl")
using .ModiaMathTest


#--------------------------------------- Generic (model independent) definitions

mutable struct Variable
    typ
    start
    fixed::Union{Bool,Nothing}
    info::AbstractString

    Variable(;typ=Real, start=nothing, fixed=nothing, info="") = new(typ, start, fixed, info)
end
const VariableDict = Dict{Symbol, Variable}


mutable struct ModelResult
    t::AbstractVector
    model::AbstractVector
    variables::VariableDict

    ModelResult(m) = new(typeof(m.simulationState.time)[],
                         typeof(m.model)[],
                         m.variables)
end

function addToResult!(result::ModelResult, t, m)::Nothing
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


mutable struct PendulumGeneric{StateType,RealType}
    # Required:
    simulationState::SimulationState{StateType, RealType}

    # Defined by model
    model::Model{StateType,RealType}   # Data structure of model variable values
    variables::VariableDict            # Data structure of model variable definition
end
const Pendulum = PendulumGeneric{Float64, Float64}


# Only needed if the same problem is used several times in solve(..)
# However in ModiaMath, problem is always newly constructed for solve.
#
# The following definition tells DifferentialEquations.jl that a SimulationModel is mutable (otherwise stack overflow error for two calls of solve!)
# using ArrayInterface
#ArrayInterface.ismutable(::Type{PendulumGeneric}) = true
#ArrayInterface.ismutable(::PendulumGeneric) = true



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
@time solution = solve!(pendulum, DifferentialEquations.Tsit5(), stopTime=7.0, tolerance=1e-4)
standardPlot(solution, 1, 2)

# DAE integrator with FLoat64
      solution = solve!(pendulum, DifferentialEquations.IDA(), stopTime=7.0, tolerance=1e-4)
@time solution = solve!(pendulum, DifferentialEquations.IDA(), stopTime=7.0, tolerance=1e-4)
standardPlot(solution, 3, 4)

# ODE integrator with DoubleFloats
pendulum = PendulumGeneric{DoubleFloat{Float64},Float64}(m=1.1)
      solution = solve!(pendulum, DifferentialEquations.Tsit5(), stopTime=7.0, tolerance=1e-4)
@time solution = solve!(pendulum, DifferentialEquations.Tsit5(), stopTime=7.0, tolerance=1e-4)
standardPlot(solution, 5, 6)

# ODE integrator with Measurements
pendulum = PendulumGeneric{Measurement{Float64},Float64}(L=1.0±0.1, m=1.0±0.1)
      solution = solve!(pendulum, DifferentialEquations.Tsit5(), stopTime=7.0, tolerance=1e-4)
@time solution = solve!(pendulum, DifferentialEquations.Tsit5(), stopTime=7.0, tolerance=1e-4)
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