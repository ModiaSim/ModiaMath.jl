# License for this file: MIT (expat)
# Copyright 2020, DLR Institute of System Dynamics and Control

"""
    module ModiaMathTest_ode_pendulum2

ODE-model of a mass point attached via a rod to a revolute joint with 2 states.

This test is used to evalute how to use DifferentialEquations.jl in ModiaMath.
The difference to ModiaMathTest_ode_pendulum.jl is that the model just uses a Float64 type.
"""
module ModiaMathTest_ode_pendulum2

import DifferentialEquations
using PyPlot        # For plotting

using StaticArrays  # Test of StaticArrays with parameterized data


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
    t::Vector{Float64}
    model::AbstractVector
    variables::VariableDict

    ModelResult(m) = new(Float64[], typeof(m.model)[], m.variables)
end


function addToResult!(result::ModelResult, t, m)::Nothing
    push!(result.t, deepcopy(t))
    push!(result.model, deepcopy(m.model))
    return nothing;
end


#---------------------------------------- Specific to a particular model (here: Pendulum)

mutable struct Model
    # Parameters
    L::Float64
    m::Float64
    d::Float64
    g::Float64

    # Variables during integration
    phi::Float64
    w::Float64
    der_w::Float64

    # Variables at communication points
    r::SVector{2,Float64}   # absolute position of end point

    function Model(;L=1.0, m=1.0, d=0.1, g=9.81)
        @assert(L > 0)
        @assert(m > 0)
        @assert(d >= 0)
        new(L, m, d, g, 0, 0, 0, SVector{2,Float64}(0,0))
    end
end


mutable struct Pendulum
    # Required:
    simulationState::SimulationState{Float64,Float64}

    # Defined by model
    model::Model               # Data structure of model variable values
    variables::VariableDict    # Data structure of model variable definition
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


function InitialStates(model::Model, var::VariableDict)::Vector{Float64}
    x_start = [Float64(var[:phi].start),
               Float64(var[:w].start)]
    return x_start
end

function setStates!(obj::Pendulum, x)::Nothing
    m = obj.model
    m.phi = x[1]
    m.w   = x[2]
    return nothing
end

function getDerivatives!(obj::Pendulum, der_x)::Nothing
    m = obj.model
    der_x[1] = m.w
    der_x[2] = m.der_w
    return nothing
end

function evaluate!(obj::Pendulum)::Nothing
    m   = obj.model
    sim = obj.simulationState
    m.der_w = (-m.m * m.g * m.L * sin(m.phi) - m.d * m.w) / (m.m * m.L^2)

    if sim.storeResult
        m.r = SVector{2,Float64}(m.L*sin(m.phi), -m.L*cos(m.phi))
    end
    return nothing
end


"Pendulum: Constructor for Pendulum SimulationModel"
function Pendulum(;kwargs...)
    model     = Model(;kwargs...)
    variables = ModelVariables()
    x_start   = InitialStates(model, variables)

    simulationState = SimulationState{Float64,Float64}("Pendulum", x_start, setStates!, getDerivatives!, evaluate!, ModelResult, addToResult!)
    simulationModel = Pendulum(simulationState, model, variables)
end




#------------------------------------------------------------------------------------------------
#---------------------------------------- Simulation of the pendulum model
#------------------------------------------------------------------------------------------------

function standardPlot(solution, fig1, fig2)
    result = solution.prob.p.simulationState.result
    Float64 = string( typeof( solution.prob.p.simulationState.x_start[1] ) )
    figure(fig1)
    clf()
    t   = result.t
    phi = getfield.(result.model, :phi)
    plot(t, phi, label="\$\\varphi\$")
    grid(true)
    legend()
    title("SimplePendulum with " * Float64)

    figure(fig2)
    clf()
    r  = getfield.(result.model, :r)
    r1 = getindex.(r,1)
    r2 = getindex.(r,2)
    plot(t, r1, label="r[1]")
    plot(t, r2, label="r[2]")
    grid(true)
    legend()
    title("SimplePendulum with " * Float64)
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


end