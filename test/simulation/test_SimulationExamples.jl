# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

module test_SimulationExamples

import ModiaMath

# Desired:
#   using Test
#
# In order that Test needs not to be defined in the user environment, it is included via ModiaMath:
using ModiaMath.Test


include(joinpath(ModiaMath.path, "examples", "Simulate_Pendulum.jl"))
import .Simulate_Pendulum
pendulumModel = Simulate_Pendulum.Pendulum(L=0.5, m=1.0, d=0.1, phi0_deg=60.0)
pendulum1 = ModiaMath.SimulationModel(pendulumModel; structureOfDAE=ModiaMath.DAE_NoSpecialStructure)
pendulum2 = ModiaMath.SimulationModel(pendulumModel; structureOfDAE=ModiaMath.DAE_ExplicitDerivatives)
pendulum3 = ModiaMath.SimulationModel(pendulumModel; structureOfDAE=ModiaMath.DAE_LinearDerivativesAndConstraints)

@testset "ModiaMath: examples/Simulate_Pendulum.jl" begin 
    result = ModiaMath.simulate!(pendulum1, stopTime=10.0, tolerance=1e-8, interval=0.1, log=true) 
    phi = result.series["phi"]
    w   = result.series["w"]

    @test isapprox(phi[end], 0.120289; atol=1e-3 )
    @test isapprox(w[end], 0.273965; atol=1e-2 )


    result = ModiaMath.simulate!(pendulum2, stopTime=10.0, tolerance=1e-8, interval=0.1, log=true) 
    phi = result.series["phi"]
    w   = result.series["w"]

    @test isapprox(phi[end], 0.120289; atol=1e-3 )
    @test isapprox(w[end], 0.273965; atol=1e-2 )


    result = ModiaMath.simulate!(pendulum3, stopTime=10.0, tolerance=1e-8, interval=0.1, log=true) 
    phi = result.series["phi"]
    w   = result.series["w"]

    @test isapprox(phi[end], 0.120289; atol=1e-3 )
    @test isapprox(w[end], 0.273965; atol=1e-2 )
end


include(joinpath(ModiaMath.path, "examples", "Simulate_FreeBodyRotation.jl"))
import .Simulate_FreeBodyRotation.FreeBodyRotation
import .Simulate_FreeBodyRotation.result

@testset "ModiaMath: examples/Simulate_FreeBodyRotation.jl" begin 
    simulationModel = ModiaMath.SimulationModel(FreeBodyRotation(), stopTime=5.0, tolerance=1e-8)
    result2         = ModiaMath.simulate!(simulationModel, log=true)

    q1 = result.series["q"]    
    w1 = result.series["w"]
    q2 = result2.series["q"]    
    w2 = result2.series["w"]

    @test isapprox(q1[end,:], q2[end,:]; atol=1e-3)
    @test isapprox(w1[end,:], w2[end,:]; atol=1e-3)    
end


end