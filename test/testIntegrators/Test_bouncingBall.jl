module Test_bouncingBall

#=
ContinuousCallback: affect! is not always called at the crossed condition

https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/551


Event restart: Stepsize not adapted for Tsit5 #558

https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/558
=#


using DifferentialEquations, Sundials, PyPlot

function f(du,u,p,t)
    du[1] = u[2]
    du[2] = -9.81
end

function condition(u,t,integrator) # Event when event_f(u,t) == 0
    z = u[1] + 1e-12
    # println("condition called at time = ", t, ", z = ", z)
    return z
end

function affect_neg(integrator)
    println("Event (affect_neg) at time = ", integrator.t, ", h = ", integrator.u[1])
    integrator.u[2] = -0.8*integrator.u[2]
    #set_proposed_dt!(integrator, 1e-5)
    #println("    after set_proposed_dt!: integrator.dt = ", integrator.dt)
    println("   before reset: integrator.dt = ", integrator.dt)
    auto_dt_reset!(integrator)
    println("   after reset : integrator.dt = ", integrator.dt)
    set_proposed_dt!(integrator, integrator.dt)
end

cb = ContinuousCallback(condition,nothing,affect_neg! = affect_neg)

u0 = [1.0,0.0]
tspan = (0.0,3.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5(),saveat=0.01,callback=cb)
#sol = solve(prob,CVODE_BDF(),saveat=0.01,callback=cb)
#=
println("propertynames(prob) = ", propertynames(prob))
println("propertynames(sol) = ", propertynames(sol))
println("propertynames(sol.prob) = ", propertynames(sol.prob))
println("typeof(sol.prob.problem_type) = ", typeof(sol.prob.problem_type))
println("typeof(sol.alg) = ", typeof(sol.alg))
=#

figure(1)
clf()
plot(sol.t, getindex.(sol.u, 1), label="h")

end