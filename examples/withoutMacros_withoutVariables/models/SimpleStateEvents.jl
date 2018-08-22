# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module SimpleStateEvents

Model of a mass controlled by a two-point controller leading to state events.
  
# Modelica model

    f = if s > 0 then 0 else fMax;
    v = der(s);
    m*der(v) + d*v + k*s = f 
"""
module SimpleStateEvents

import ModiaMath


mutable struct Model <: ModiaMath.AbstractSimulationModel
    simulationState::ModiaMath.SimulationState
   
    # Parameters   
    m::Float64
    k::Float64
    d::Float64
    fMax::Float64
      
    function Model(;m=1.0, k=1.0, d=0.1, fMax=1.5, s0=2.0)
        @assert(m > 0.0)
        @assert(k > 0.0)
        @assert(d >= 0.0)
        @assert(fMax > 0.0)
        simulationState = ModiaMath.SimulationState("SimpleStateEvents", getModelResidues!, [s0;0.0], getVariableName;
                                nz=1, nw=1)            
        new(simulationState, m, k, d, fMax)
    end
end 

getVariableName(model, vcat, vindex) = ModiaMath.getVariableName(model, vcat, vindex;
                                                         xNames=["s", "v"],
                                                         wNames=["f"])
   
function getModelResidues!(m::Model, t::Float64, _x::Vector{Float64}, _derx::Vector{Float64}, _r::Vector{Float64}, _w::Vector{Float64})  
    sim = m.simulationState
  
    s = _x[1]
    v = _x[2]
   
    sPos = ModiaMath.positive!(sim, 1, s, "s")
    f    = sPos ? 0.0 : m.fMax 
    derv = (f - m.d * v - m.k * s) / m.m
   
    _r[1] = v    - _derx[1]
    _r[2] = derv - _derx[2]
   
    if ModiaMath.isStoreResult(m.simulationState)
        _w[1] = f
    end
    return nothing
end
 
end