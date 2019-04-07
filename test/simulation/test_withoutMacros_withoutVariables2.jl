# License for this file: MIT (expat)
# Copyright 2016-2018, DLR Institute of System Dynamics and Control


"""
    module test_withoutMacros_withoutVariables2

Test function ModiaMath.simulate(..) with models that have events
"""
module  test_withoutMacros_withoutVariables2

import ModiaMath

# Desired:
#   using Test
#
# In order that Test needs not to be defined in the user environment, it is included via ModiaMath:
using ModiaMath.Test



include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "models", "SimpleStateEvents.jl"))
import .SimpleStateEvents

include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "models", "BouncingBall.jl"))
import .BouncingBall

include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "models", "IdealClutch.jl"))
import .IdealClutch


@testset "\nTest ModiaMath: withoutMacros_withoutVariables/*.jl with events" begin

    @testset "Simulate SimpleStateEvents" begin
        m1 = SimpleStateEvents.Model()
        result = ModiaMath.simulate!(m1, stopTime=10.0, interval=0.1)
        time = result["time"]
        f    = result["f"]
        nStateEvents = m1.simulationState.statistics.nStateEvents
        tEvents = [1.62283; 3.36771; 6.51321; 8.04188]
        iEvents = [16 + 2, 33 + 2 + 2, 65 + 2 + 4, 80 + 2 + 6]

        @test length(time) == 101 + 4 * 2
        @test nStateEvents == 4
        @test isapprox(tEvents, time[iEvents]; atol=0.001 )
    end

    @testset "Simulate BouncingBall" begin
        m2 = BouncingBall.Model()
        result = ModiaMath.simulate!(m2, stopTime=3.0, interval=0.01)
        nStateEvents = m2.simulationState.statistics.nStateEvents

        @test nStateEvents == 18
    end

    @testset "Simulate IdealClutch" begin
        m2 = IdealClutch.Model()
        result = ModiaMath.simulate!(m2, stopTime=500.0)

        w1_end = result["inertia1.w"][end]
		w2_end = result["inertia2.w"][end]
		w_end_required = 38.9277466565
		
        @test isapprox(w1_end, w_end_required; atol=0.001 )
        @test isapprox(w2_end, w_end_required; atol=0.001 )		 
    end
	
end

include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "Simulate_SimpleStateEvents.jl"))
include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "Simulate_BouncingBall.jl"))
include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "Simulate_IdealClutch.jl"))

end