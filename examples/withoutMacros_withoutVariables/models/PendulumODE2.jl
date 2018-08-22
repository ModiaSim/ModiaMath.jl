# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module PendulumODE2

ODE-model of a mass point attached via a rod to a revolute joint with 2 states.
"""
module PendulumODE2

import ModiaMath

     
mutable struct Model <: ModiaMath.AbstractSimulationModel
    simulationState::ModiaMath.SimulationState
   
    # Parameters
    L::Float64
    m::Float64
    d::Float64
    g::Float64

    function Model(;L=1.0, m=1.0, d=0.1, g=9.81, phi0_deg=90)
        @assert(L > 0.0)
        @assert(m > 0.0)
        @assert(d >= 0.0)
        simulationState = ModiaMath.SimulationState("PendulumODE2", getModelResidues!, [phi0_deg * pi / 180; 0.0], getVariableName; nw=2)
        new(simulationState, L, m, d, g)
    end
end 

getVariableName(model, vcat, vindex) = ModiaMath.getVariableName(model, vcat, vindex;
                                                         xNames=["phi"    , "w"],
                                                         derxNames=["der_phi", "der_w"],
                                                         wNames=["rx", "ry"]) 

function getModelResidues!(m::Model, t::Float64, _x::Vector{Float64}, _derx::Vector{Float64}, _r::Vector{Float64}, _w::Vector{Float64})  
    phi  = _x[1]
    w    = _x[2]  
    derw = _derx[2]
    _r[1] = w    - _derx[1]
    _r[2] = m.m * m.L^2 * derw + (m.m * m.g * m.L * sin(phi) + m.d * w)
   
    if ModiaMath.isStoreResult(m.simulationState)
        rx =  m.L * sin(phi)
        ry = -m.L * cos(phi)
        _w[1] = rx
        _w[2] = ry
    end
    return nothing
end

end