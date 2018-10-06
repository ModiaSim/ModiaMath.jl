# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

module test_SimulationExamples

import ModiaMath

# Desired:
#   using Test
#
# In order that Test needs not to be defined in the user environment, it is included via ModiaMath:
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using ModiaMath.Test
end


include(joinpath(ModiaMath.path, "examples", "Simulate_Pendulum.jl"))
import .Simulate_Pendulum
pendulum = ModiaMath.SimulationModel(Simulate_Pendulum.Pendulum(L=0.5, m=1.0, d=0.1, phi0_deg=60.0))

@testset "ModiaMath: examples/Simulate_Pendulum.jl" begin 
    result = ModiaMath.simulate!(pendulum, stopTime=10.0, tolerance=1e-8, interval=0.1, log=true) 
    phi = result.series["phi"]
    w   = result.series["w"]

    @test isapprox(phi[end], 0.120289; atol=1e-3 )
    @test isapprox(w[end], 0.273965; atol=1e-2 )
end


include(joinpath(ModiaMath.path, "examples", "Simulate_FreeBodyRotation.jl"))


end