# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module Simulate_FreeBodyRotation  (examples/withoutMacros)

Simulate a free body rotation with quaternions (version that uses Variables, but no macros)

A body can rotate freely around the origin and a sin-wave torque vector is applied on it.
The system is mathematically described as a fully-implicit, index-1 DAE model of a body
with unconstrained rotational motion in 3D using quaternions.

A simulation tool of a corresponding Modelica model that transforms the model to an ODE
would need to use dynamic state selection and solve nonlinear algebraic equations
within the ODE. Due to the implicit index-1 DAE description, this is much simpler here.

# DAE equations
```
  t := A.*sin(2*pi*freqHz*time + phase);
  0  = w - 2*([q[4],  q[3], -q[2], -q[1];
              -q[3],  q[4],  q[1], -q[2];
             q[2], -q[1],  q[4], -q[3]])*der(q);
  0  = I*der(w) + cross(w, I*w) - t;
  0  = dot(q,q)-1.0;

Initialization:
  q = [0.1, 0.5, 0.0, 1.0]   # inconsistent values (fixed during initialization)
  w = zeros(3)

Reference Modelica model:
   ModelicaReferenceModels.ODAEs.FreeBody
```
"""
module Simulate_FreeBodyRotation

using ModiaMath
using ModiaMath.LinearAlgebra  # included via ModiaMath, to avoid requirement to add it in the standard environment
using ModiaMath.StaticArrays   # included via ModiaMath, to avoid requirement to add it in the standard environment


#            q[1] = 0.1 changed to 0.08908708063747484
#            q[2] = 0.5 changed to 0.445435403187374
#            q[4] = 1.0 changed to 0.8908708063747479
#                             q0     = [0.1, 0.5, 0.0, 1.0],
# q0     = [0.08908708063747484, 0.445435403187374, 0.0, 0.8908708063747479],

mutable struct FreeBodyRotationWithoutMacro <: ModiaMath.AbstractComponentWithVariables
    _internal::ModiaMath.ComponentInternal

    # Parameters
    I::SMatrix{3,3,Float64,9}
    A::SVector{3,Float64}
    freqHz::SVector{3,Float64}
    phase::SVector{3,Float64}
    q0::SVector{4,Float64}
    w0::SVector{3,Float64}

    # Variables
    q::RealSVector{4}
    derq::RealSVector{4}
    w::RealSVector3
    derw::RealSVector3
    tau::RealSVector3
    residue_w::RealSVector3
    residue_t::RealSVector3
    residue_q::RealScalar

    function FreeBodyRotationWithoutMacro(;I=Diagonal([1.0,2.0,3.0]),
                                         A=[3.0,4.0,5.0],
                                         freqHz=[0.3,0.2,0.1],
                                         phase=[0,0.5235987755983,1.0471975511966],
                                         q0=[0.1, 0.5, 0.0, 1.0],
                                         w0=zeros(3))
        this = new(ModiaMath.ComponentInternal(:FreeBodyRotationWithoutMacro, nothing), I, A, freqHz, phase, q0, w0)

        @assert(size(I) == (3, 3))
        @assert(isposdef(I))
        @assert(length(A) == 3)
        @assert(length(freqHz) == 3)
        @assert(length(phase) == 3)
        @assert(minimum(A) > 0.0)   
    
        q         = RealSVector{4}(:q, this, numericType=ModiaMath.XD_IMP, info="Quaternions",                         start=q0, fixed=false)
        derq      = RealSVector{4}(:derq, this, numericType=ModiaMath.DER_XD_IMP, integral=q, info="der(q)",          unit="1/s")
        w         = RealSVector3(:w, this, numericType=ModiaMath.XD_IMP, info="Angular velocity",     unit="rad/s",  start=w0, fixed=true)
        derw      = RealSVector3(:derw, this, numericType=ModiaMath.DER_XD_IMP, integral=w, info="Angular acceleration", unit="rad/s^2")
        tau       = RealSVector3(:tau, this, numericType=ModiaMath.WR, info="Driving torque",       unit="N*m")
     
        residue_w = RealSVector3(:residue_w, this, numericType=ModiaMath.FD_IMP, info="Angular velocity residue")
        residue_t = RealSVector3(:residue_t, this, numericType=ModiaMath.FD_IMP, info="Angular momentum equation residue")
        residue_q = RealScalar(:residue_q, this, numericType=ModiaMath.FC, info="Quaternion constraint residue")

        return this
    end
end

function ModiaMath.computeVariables!(b::FreeBodyRotationWithoutMacro, sim::ModiaMath.SimulationState)
    time::Float64 = ModiaMath.getTime(sim)
    I             = b.I
    A             = b.A
    freqHz        = b.freqHz
    phase         = b.phase
    q             = b.q.value
    derq          = b.derq.value
    w             = b.w.value
    derw          = b.derw.value

    b.tau.value       = A .* sin.(2.0 * pi * freqHz * time + phase)
    b.residue_w.value = w - 2.0 * ( @SMatrix( [ q[4]  q[3] -q[2] -q[1];
                                            -q[3]  q[4]  q[1] -q[2];
                                             q[2] -q[1]  q[4] -q[3]] ) * derq)
    b.residue_t.value = I * derw + cross(w, I * w) - b.tau.value
    b.residue_q.value = dot(q, q) - 1.0

    return nothing
end

simulationModel = ModiaMath.SimulationModel(FreeBodyRotationWithoutMacro(), stopTime=5.0, tolerance=1e-8)
#ModiaMath.print_ModelVariables(simulationModel)
result          = ModiaMath.simulate!(simulationModel, log=true)
ModiaMath.plot(result, [:q, :w, :derw, :tau])

end
