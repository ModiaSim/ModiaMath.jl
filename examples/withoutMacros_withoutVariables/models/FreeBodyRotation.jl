# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module FreeBodyRotation

DAE model of a body with unconstrained rotational motion in 3D (using quaternions;)
    
DAE equations:
  u := A.*sin(2*pi*freqHz*time + phase);
  0 = w - 2*([Q[4],  Q[3], -Q[2], -Q[1];
             -Q[3],  Q[4],  Q[1], -Q[2];
              Q[2], -Q[1],  Q[4], -Q[3]])*der(Q);
  0 = I*der(w) + cross(w, I*w) - u;
  0 = dot(Q,Q)-1.0;

Initialization (inconsistent values)
  Q[1] = 0.1
  Q[2] = 0.5
  Q[3] = 0
  Q[4] = 1.0
  w = zeros(3)  
   
Arguments of getModelResidues! function:
   x    = [Q     ; w]
   derx = [der(Q); der(w)]
   
For testing:
   Modelica model ModelicaReferenceModels.ODAEs.FreeBody
"""
module FreeBodyRotation

import ModiaMath
using  ModiaMath.LinearAlgebra    # included via ModiaMath, to avoid requirement to add it in the standard environment


# Desired:
#   using StaticArrays
#
# In order that StaticArrays need not to be defined in the user environment, it is included via ModiaMath:
using ModiaMath.StaticArrays


mutable struct Model <: ModiaMath.AbstractSimulationModel
    simulationState::ModiaMath.SimulationState
   
    # Parameters
    m::Float64
    I::SMatrix{3,3,Float64,9}
    A::SVector{3,Float64}
    freqHz::SVector{3,Float64}
    phase::SVector{3,Float64}

    function Model(;m=1.0, I=Diagonal([1.0,2.0,3.0]), A=[3.0,4.0,5.0], freqHz=[0.3,0.2,0.1], 
                  phase=[0,0.5235987755983,1.0471975511966],
                  Q0=[0.1, 0.5, 0.0, 1.0],
                  w0=zeros(3), linearDerivatives=true)
        @assert(m > 0.0)
        @assert(size(I) == (3, 3))
        @assert(isposdef(I))
        @assert(length(A) == 3)
        @assert(length(freqHz) == 3)
        @assert(length(phase) == 3)
        @assert(minimum(A) > 0.0)
        @assert(length(Q0) == 4)
        @assert(length(w0) == 3)
        is_constraint    = fill(false,7)
        is_constraint[7] = true
        simulationState = ModiaMath.SimulationState("FreeBodyRotation", getModelResidues!, Vector{Float64}([Q0;w0]), getVariableName; 
                                                    is_constraint = is_constraint,  
                                                    structureOfDAE = linearDerivatives ? ModiaMath.DAE_LinearDerivativesAndConstraints :
                                                                                         ModiaMath.DAE_NoSpecialStructure) 
        new(simulationState, m, SMatrix{3,3,Float64,9}(I), SVector{3,Float64}(A),
                                                           SVector{3,Float64}(freqHz),
                                                           SVector{3,Float64}(phase))
    end
end 

getVariableName(model, vcat, vindex) = ModiaMath.getVariableName(model, vcat, vindex;
                                                               xNames=["Q[1]", "Q[2]", "Q[3]", "Q[4]", "w[1]", "w[2]", "w[3]"])
                                                   
function getModelResidues!(m::Model, _t::Float64, _x::Vector{Float64}, _derx::Vector{Float64}, _r::Vector{Float64}, _w::Vector{Float64})
    u = m.A .* sin.(2.0 * pi * m.freqHz * _t + m.phase)   
 
    Q    = SVector{4,Float64}(_x[1:4])
    w    = SVector{3,Float64}(_x[5:7])
    derQ = SVector{4,Float64}(_derx[1:4])
    derw = SVector{3,Float64}(_derx[5:7])

    _r[1:3] = w - 2.0 * ([ Q[4]  Q[3] -Q[2] -Q[1]; 
                       -Q[3]  Q[4]  Q[1] -Q[2];
                        Q[2] -Q[1]  Q[4] -Q[3]] * derQ)
    _r[4:6] = m.I * derw + cross(w, m.I * w) - u
    _r[7]   = dot(Q, Q) - 1.0
 
    return nothing
end

end