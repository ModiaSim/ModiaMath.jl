module Test_bouncingBall_with_measurement

#=
Strange result of multiplication

https://github.com/JuliaPhysics/Measurements.jl/issues/64
=#

using DifferentialEquations, Measurements, PyPlot

function f(du,u,p,t)
    du[1] = u[2]
    du[2] = -9.81
end

function condition(u,t,integrator)
    z = u[1] + 1e-12
    return z
end

const e = 0.7 ± 0.1

function affect_neg(integrator)
    println("Event at time = ", integrator.t, ", h = ", integrator.u[1])

    v_before = integrator.u[2]
    println("   v_before = ", v_before, ", e = ", e)
    v_after = -e*v_before
    println("   v_after  = ", v_after)

    integrator.u[2] = v_after
    auto_dt_reset!(integrator)
    set_proposed_dt!(integrator, integrator.dt)
end

cb = ContinuousCallback(condition,nothing,affect_neg! = affect_neg)

u0 = [1.0 ± 0.2, 0.0]
tspan = (0.0,2.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5(),saveat=0.01,callback=cb)

println("\n... Separate computation")
v_before = -3.1 ± 0.44
println("   v_before = ", v_before, ", e = ", e)
v_after  = -e*v_before
println("   v_after  = ", v_after)


figure(1)
clf()
plot(sol.t, getindex.(sol.u, 1), label="h")

end