module Test_bouncingBall_no_additional_values_at_events

#=
SavedValues not stored at events #557

https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/557
=#



using DifferentialEquations, PyPlot

function f(du,u,p,t)
  du[1] = u[2]
  du[2] = -9.81
end

function outputs(u, t, integrator) # called at communication points
    return 0.5*u[1]
end

function condition(u,t,integrator) # Event when event_f(u,t) == 0
 z = u[1] + 1e-12
 return z
end

function affect_neg(integrator)
  integrator.u[2] = -0.8*integrator.u[2]
end

u0 = [1.0,0.0]
tspan = (0.0,3.0)

t_start = 0.0
t_inc   = 0.01
t_end   = 3.0
tspan   = (t_start, t_end)
tspan2  = t_start:t_inc:t_end

saved_values=SavedValues(Float64, Float64)
cb1 = SavingCallback(outputs, saved_values, saveat=tspan2)
cb2 = ContinuousCallback(condition,nothing,affect_neg! = affect_neg)

prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5(),saveat=t_inc,callback=CallbackSet(cb1,cb2))

figure(1)
clf()
plot(sol.t, getindex.(sol.u, 1), label="h")
plot(saved_values.t, saved_values.saveval, label="outputs")
grid(true)
legend()

println("length(sol.t)          = ", length(sol.t))
println("length(saved_values.t) = ", length(saved_values.t))
end