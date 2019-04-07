# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module PT1

ODE-model of PT1 block (T*dx/dt + x = 0).
"""
module PT1

import ModiaMath

   
mutable struct Model <: ModiaMath.AbstractSimulationModel
    simulationState::ModiaMath.SimulationState
   
    # Parameters
    T::Float64
   
    # Initial conditions
    x0::Float64
    
    function Model(;T=1.0, x0=1.0)
        @assert(T > 0.0)
        simulationState = ModiaMath.SimulationState("PT1", getModelResidues!, [x0];
                                                    structureOfDAE = ModiaMath.DAE_ExplicitDerivatives)
        new(simulationState, T, x0)
    end
end 
   
function getModelResidues!(m::Model, t::Float64, _x::Vector{Float64}, _derx::Vector{Float64}, _r::Vector{Float64}, _w::Vector{Float64})  
    x     = _x[1]
    derx  = -x / m.T
    _r[1] = derx - _derx[1]
    return nothing
end


end