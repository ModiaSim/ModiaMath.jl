# License for this file: MIT (expat)
# Copyright 2016-2018, DLR Institute of System Dynamics and Control

"""
    module test_withoutMacros_withoutVariables1

Test function ModiaMath.simulate(..) with standard models (no events, no sparse Jacobian),
and without Variables and without Macros
"""
module test_withoutMacros_withoutVariables1

import ModiaMath

# Desired:
#   using Test
#
# In order that Test needs not to be defined in the user environment, it is included via ModiaMath:
using ModiaMath.Test


include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "models", "PT1.jl"))
import .PT1

include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "models", "PendulumODE.jl"))
import .PendulumODE

include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "models", "PendulumDAE.jl"))
import .PendulumDAE

include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "models", "FreeBodyRotation.jl"))
import .FreeBodyRotation

include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "models", "StateSelection.jl"))
import .StateSelection

@testset "ModiaMath: withoutMacros_withoutVariables/*.jl without events" begin
    # Simulate PT1 block "T*der(x) + x = 0" in different variants
    m = PT1.Model(T=2.0, x0=1.5)
    interval = 0.1
    stopTime = 10.0
    t = collect(0.0:interval:stopTime)
    nt = length(t)

    # Analytic solution
    function fExact(m::PT1.Model, t::Vector{Float64})
        nt = length(t)
        xExact     = zeros(nt)
        der_xExact = zeros(nt)
        for i = 1:nt
            xExact[i,1]     = exp(-t[i] / m.T) * m.x0
            der_xExact[i,1] = -xExact[i] / m.T
        end
        return (xExact, der_xExact)
    end

    # Test whether simulation results in the exact solution
    @testset "Simulate PT1 with default tolerance" begin
        tolerance = 1e-4
        result    = ModiaMath.simulate!(m, stopTime=stopTime, interval=interval, tolerance=tolerance)
        (xExact, der_xExact) = fExact(m, t)

        @test isapprox(result["x[1]"], xExact    ; atol=10 * tolerance )
        @test isapprox(result["der(x[1])"], der_xExact; atol=100 * tolerance )
        @test m.simulationState.statistics.sparseSolver == false 

    end
   
    @testset "Simulate PT1 with tolerance = 1e-8 and log=true" begin
        tolerance = 1.0e-8
        result = ModiaMath.simulate!(m, stopTime=stopTime, log=true, interval=interval, tolerance=tolerance)
        (xExact, der_xExact) = fExact(m, t)

        @test isapprox(result["x[1]"], xExact    ; atol=10 * tolerance )
        @test isapprox(result["der(x[1])"], der_xExact; atol=100 * tolerance )
        @test m.simulationState.statistics.sparseSolver == false
    end

    @testset "Simulate PendulumODE" begin
        tolerance = 1e-8
        m3 = PendulumODE.Model(L=0.5, m=1.0, d=0.1, phi0_deg=60.0)   # L,m,d,phi_deg
        result = ModiaMath.simulate!(m3, stopTime=stopTime, interval=interval, tolerance=tolerance)
        phi = result["phi"]
        w   = result["w"]

        @test isapprox(phi[end], 0.120289; atol=1e-3 )
        @test isapprox(w[end], 0.273965; atol=1e-2 )
    end

    @testset "Simulate PendulumDAE" begin
        stopTime  = 2.0
        tolerance = 1e-8
        m4 = PendulumDAE.Model()
        result = ModiaMath.simulate!(m4, stopTime=stopTime, interval=interval, tolerance=tolerance)
        x  = result["x"]
        y  = result["y"]
        vx = result["vx"]
        vy = result["vy"]
        dervx = result["der(vx)"]
        dervy = result["der(vy)"]

        @test isapprox(x[end],  0.542797; atol=1e-5 )
        @test isapprox(y[end], -0.839864; atol=1e-5 )
        @test isapprox(vx[end],  1.80223 ; atol=1e-5 )
        @test isapprox(vy[end],  1.16476 ; atol=1e-5 )
        @test isapprox(dervx[end], -6.97155 ; atol=1e-5 )
        @test isapprox(dervy[end],  0.977022; atol=1e-5 )
    end
   

    @testset "Simulate PendulumDAE with x_fixed=true" begin
        stopTime  = 2.0
        tolerance = 1e-8
        m4 = PendulumDAE.Model(x_fixed=true)
        result = ModiaMath.simulate!(m4, stopTime=stopTime, interval=interval, tolerance=tolerance)
        x  = result["x"]
        y  = result["y"]
        vx = result["vx"]
        vy = result["vy"]
        dervx = result["der(vx)"]
        dervy = result["der(vy)"]

        @test isapprox(x[end],  0.431501; atol=1e-5 )
        @test isapprox(y[end], -0.902112; atol=1e-5 )
        @test isapprox(vx[end],  1.2889  ; atol=1e-5 )
        @test isapprox(vy[end],  0.616512; atol=1e-5 )
        @test isapprox(dervx[end], -4.69952 ; atol=1e-5 )
        @test isapprox(dervy[end],  0.01497 ; atol=1e-5 )
    end
   
    @testset "Simulate FreeBodyRotation" begin
        stopTime  = 5.0
        tolerance = 1e-8
        m5 = FreeBodyRotation.Model()
        result = ModiaMath.simulate!(m5, stopTime=stopTime, interval=interval, tolerance=tolerance)

        #=
        @test isapprox(x[end]    ,  0.542797; atol=1e-3 )
        @test isapprox(y[end]    , -0.839864; atol=1e-3 )
        @test isapprox(vx[end]   ,  1.80223 ; atol=1e-3 )
        @test isapprox(vy[end]   ,  1.16476 ; atol=1e-3 )
        @test isapprox(dervx[end], -6.97155 ; atol=1e-3 )
        @test isapprox(dervy[end],  0.97702 ; atol=1e-3 )
        =#
    end

    @testset "Simulate State Selection" begin
        m6 = StateSelection.Model()
        result = ModiaMath.simulate!(m6, stopTime=1.0)
		f3 = result["f[3]"][end]

        @test isapprox(f3[end], 9.81; atol=1e-3 )
    end   
end 

include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "Simulate_PT1.jl"))
include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "Simulate_PendulumODE.jl"))
include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "Simulate_PendulumDAE.jl"))
include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "Simulate_FreeBodyRotation.jl"))
include(joinpath(ModiaMath.path, "examples", "withoutMacros_withoutVariables", "Simulate_StateSelection.jl"))

end