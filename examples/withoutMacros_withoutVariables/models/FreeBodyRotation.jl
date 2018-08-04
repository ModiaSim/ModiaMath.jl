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

     
mutable struct Model <: ModiaMath.AbstractSimulationModel
   simulationState::ModiaMath.SimulationState
   
   # Parameters
   m::Float64
   I::Matrix{Float64}
   A::Vector{Float64}
   freqHz::Vector{Float64}
   phase::Vector{Float64}

   function Model(;m=1.0, I=diagm([1.0,2.0,3.0]), A=[3.0,4.0,5.0], freqHz=[0.3,0.2,0.1], 
                  phase=[0,0.5235987755983,1.0471975511966],
                  Q0::Vector{Float64}=[0.1, 0.5, 0.0, 1.0],
                  w0::Vector{Float64}=zeros(3))
      @assert(m > 0.0)
      @assert(size(I)==(3,3))
      @assert(isposdef(I))
      @assert(length(A)==3)
      @assert(length(freqHz)==3)
      @assert(length(phase)==3)
      @assert(minimum(A) > 0.0)

      simulationState = ModiaMath.SimulationState("FreeBodyRotation", getModelResidues!, [Q0;w0], getVariableName;  nc = 1) 
      new(simulationState,m,I,A,freqHz,phase)
   end
end 

getVariableName(model,vcat,vindex) = ModiaMath.getVariableName(model,vcat,vindex;
                                                               xNames = ["Q[1]", "Q[2]", "Q[3]", "Q[4]", "w[1]", "w[2]", "w[3]"])
                                                   
function getModelResidues!(m::Model, _t::Float64, _x::Vector{Float64}, _derx::Vector{Float64}, _r::Vector{Float64}, _w::Vector{Float64})
   u = m.A .* sin.(2.0*pi*m.freqHz*_t + m.phase)   
 
   Q    = _x[1:4]
   w    = _x[5:7]
   derQ = _derx[1:4]
   derw = _derx[5:7]

   _r[1:3] = w - 2.0*([ Q[4]  Q[3] -Q[2] -Q[1]; 
                       -Q[3]  Q[4]  Q[1] -Q[2];
                        Q[2] -Q[1]  Q[4] -Q[3]]*derQ)
   _r[4:6] = m.I*derw + cross(w, m.I*w) - u
   _r[7]   = dot(Q,Q) - 1.0
 
   return nothing
end

end