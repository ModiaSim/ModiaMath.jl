# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module StateSelection

DAE-model to test manual state selection of simple multibody system.

This model is mainly used to test initialization
"""
module StateSelection

import ModiaMath
using  ModiaMath.LinearAlgebra

mutable struct Model <: ModiaMath.AbstractSimulationModel
    simulationState::ModiaMath.SimulationState
   
    # Parameters
    n::Vector{Float64}
    m::Float64
    g::Float64
    s0::Float64
    v0::Float64

    function Model(;n=[0.0, 1.0, 0.0], m=1.0, g=9.81, s0=0.0, v0=0.0)
        @assert(length(n) == 3)
        @assert(m > 0.0)
        simulationState = ModiaMath.SimulationState("StateSelection", getModelResidues!, [s0,0.0,0.0,0.0,v0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], getVariableName;
                                x_fixed=[true, false,false,false, true, false,false,false, false, false,false,false],
								structureOfDAE = ModiaMath.DAE_NoSpecialStructure,
                                hev=1e-4)
        new(simulationState, n, m, g, s0, v0)
    end
end 

   
getVariableName(model, vcat, vindex) = ModiaMath.getVariableName(model, vcat, vindex;
                                                xNames=["s", "f[1]", "f[2]", "f[3]", "sd", "der_der_r[1]","der_der_r[2]","der_der_r[3]",
                                                             "der_der_s", "der_v[1]", "der_v[2]", "der_v[3]"]) 
                   
function getModelResidues!(m::Model, t::Float64, _x::Vector{Float64}, _derx::Vector{Float64}, _r::Vector{Float64}, _w::Vector{Float64})   
    # map x
    s  = _x[1]
    f  = _x[2:4]
	g  = [0.0, 0.0, m.g]
    sd = _x[5]
    der_der_r = _x[6:8]
    der_der_s = _x[9]
    der_v     = _x[10:12]

    # map derx
    ders  = _derx[1]
    dersd = _derx[5] 
     
    # compute residues
    _r[1]    = sd - ders
    _r[2:4]  = m.n * der_der_s - der_der_r 
    _r[5:7]  = der_der_r - der_v 
    _r[8:10] = f - m.m*g - m.m * der_v
    _r[11]   = dot(m.n, f)
    _r[12]   = der_der_s - dersd
    return nothing
end

end