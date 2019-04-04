# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module BouncingBall

Model of a ball bouncing on ground.
"""
module BouncingBall

import ModiaMath

mutable struct Model <: ModiaMath.AbstractSimulationModel
    simulationState::ModiaMath.SimulationState
   
    # Parameters
    h0::Float64   # initial height
    e::Float64    # coefficient of restitution
    g::Float64    # gravity constant
   
    # Discrete variables
    flying::Bool
   
    function Model(;h0=1.0, e=0.7, g=9.81)
        @assert(h0 > 0.0)   
        @assert(0.0 <= e <= 1.0)   
        simulationState = ModiaMath.SimulationState("BouncingBall", getModelResidues!, [h0;0.0], getVariableName;
                                                    nz=1, nw=1)
        new(simulationState, h0, e, g, true)
    end
end

getVariableName(model, vcat, vindex) = ModiaMath.getVariableName(model, vcat, vindex;
                                                               xNames=["h", "v"],
                                                               wNames=["flying"])

function getModelResidues!(m::Model, t::Float64, _x::Vector{Float64}, _derx::Vector{Float64}, _r::Vector{Float64}, _w::Vector{Float64})  
    sim = m.simulationState
    if ModiaMath.isInitial(sim)
        if ModiaMath.isLogInfos(sim)
            println("        flying = ", m.flying)
        end
    end
    if ModiaMath.isTerminal(sim)
        if ModiaMath.isLogInfos(sim)
            println("\n      BouncingBall model is terminated (flying = ", m.flying, ")")
        end
        return
    end
   
    h = _x[1]
    v = _x[2]
   
    if ModiaMath.edge!(sim, 1, -h, "-h")
        v = -m.e * v
        if v < 0.01
            m.flying = false
            v = 0.0
            if ModiaMath.isLogInfos(sim)
                println("        flying = ", m.flying)
            end         
        end
        _x[2] = v    # re-initialize state vector x
    end

    derh = v
    derv = m.flying ? -m.g : 0.0

    _r[1] = derh - _derx[1]
    _r[2] = derv - _derx[2]
    _w[1] = m.flying

    return nothing
end


end