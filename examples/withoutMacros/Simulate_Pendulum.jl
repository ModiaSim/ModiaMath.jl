# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module Simulate_PendulumWithoutMacro (examples/withoutMacros)

Simulate Pendulum as ODE model (version that uses Variables, but no macros)
"""
module Simulate_Pendulum

using ModiaMath

mutable struct PendulumWithoutMacro <: ModiaMath.AbstractComponentWithVariables
    _internal::ModiaMath.ComponentInternal

    # Parameters
    L::Float64
    m::Float64
    d::Float64
    g::Float64
    phi0_deg::Float64

    # Variables
    phi::ModiaMath.RealScalar
    w::ModiaMath.RealScalar
    a::ModiaMath.RealScalar 

    function PendulumWithoutMacro(;L=1.0, m=1.0, d=0.1, g=9.81, phi0_deg=90)
        this = new(ModiaMath.ComponentInternal(:PendulumWithoutMacro, nothing))
        this.L = L
        this.m = m
        this.d = d
        this.g = g
        this.phi0_deg = phi0_deg

        @assert(L > 0.0)
        @assert(m > 0.0)
        @assert(d >= 0.0)

        phi = ModiaMath.RealScalar(start=phi0_deg * pi / 180, unit="rad", fixed=true, info="Relative rotation angle", numericType=ModiaMath.XD_EXP)
        ModiaMath.initComponent!(this, phi, :phi)

        w   = ModiaMath.RealScalar(start=0.0, unit="rad/s", fixed=true, info="Relative angular velocity", integral=phi, numericType=ModiaMath.XD_EXP)
        ModiaMath.initComponent!(this, w, :w)

        a   = ModiaMath.RealScalar(unit="rad/s^2", info="Relative angular acceleration", integral=w, numericType=ModiaMath.DER_XD_EXP) 
        ModiaMath.initComponent!(this, a, :a)

        return this
    end
end

function ModiaMath.computeVariables!(p::PendulumWithoutMacro, sim::ModiaMath.SimulationState)  
    L = p.L
    m = p.m
    d = p.d
    g = p.g
    phi = p.phi.value
    w   = p.w.value
   
    p.a.value = (-m * g * L * sin(phi) - d * w) / (m * L^2)
end

simulationModel = ModiaMath.SimulationModel(PendulumWithoutMacro(L=0.8, m=0.5, d=0.2), stopTime=5.0)
#ModiaMath.print_ModelVariables(simulationModel)
result = ModiaMath.simulate!(simulationModel, log=true) 

ModiaMath.plot(result, [:phi, :w, :a])

end