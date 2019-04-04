# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module Simulate_Pendulum

Demonstrate how a pendulum can be defined, simulated and plotted
with ModiaMath.
"""
module Simulate_Pendulum

using ModiaMath
using ModiaMath.StaticArrays   # included via ModiaMath, to avoid requirement to add it in the standard environment


@component Pendulum(;L=1.0, m=1.0, d=0.1, g=9.81, phi0_deg=90.0) begin
    @assert(L > 0.0)
    @assert(m > 0.0)
    @assert(d >= 0.0)
  
    phi = RealScalar(start=phi0_deg*pi/180, unit="rad"    , fixed=true, nominal=1.0, info="Relative rotation angle",                     numericType=ModiaMath.XD_EXP)
    w   = RealScalar(start=0.0            , unit="rad/s"  , fixed=true, nominal=1.0, info="Relative angular velocity",     integral=phi, numericType=ModiaMath.XD_EXP)
    a   = RealScalar(                       unit="rad/s^2",                          info="Relative angular acceleration", integral=w  , numericType=ModiaMath.DER_XD_EXP) 
    r   = RealSVector{2}(                   unit="m"      ,                          info="Tip position of pendulum"     ,               numericType=ModiaMath.WC)
end

function ModiaMath.computeVariables!(p::Pendulum, sim::ModiaMath.SimulationState)  
    L = p.L
    m = p.m
    d = p.d
    g = p.g
    phi = p.phi.value
    w   = p.w.value
   
    p.a.value = (-m*g*L*sin(phi) - d*w) / (m*L^2)

    if ModiaMath.isStoreResult(sim)
        p.r.value = @SVector [L*sin(phi), -L*cos(phi)]
    end
end

simulationModel = ModiaMath.SimulationModel(Pendulum(L=0.8, m=0.5, d=0.2); stopTime=5.0, structureOfDAE=ModiaMath.DAE_ExplicitDerivatives)
result = ModiaMath.simulate!(simulationModel; log=true) 

plot(result, [(:phi, :w) :a])
#plot(result, "r[2]", xAxis="r[1]", figure=2)

# ModiaMath.print_ModelVariables(simulationModel)
# println("result variables = ", result)

end