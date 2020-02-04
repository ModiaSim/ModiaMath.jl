module Test_bouncingBall_with_measurement

#=
Strange result of bouncing ball with Measurement.jl

https://github.com/JuliaPhysics/Measurements.jl/issues/64

Issue fixed in Measurements version 2.2.0
=#

using DifferentialEquations, Measurements, PyPlot


const g   = 9.81       # Gravity constant
const h0  = 1.0 ± 0.2  # Initial height
const cof = 0.7 ± 0.1  # Coefficient of restitution


function f(du,u,p,t)
    du[1] = u[2]
    du[2] = -g
end

function condition(u,t,integrator)
    z = u[1] + 1e-12
    return z
end


function affect_neg(integrator)
    println("Event at time = ", integrator.t, ", h = ", integrator.u[1])

    v_before = integrator.u[2]
    println("   v_before = ", v_before)
    v_after = -cof*v_before
    println("   v_after  = ", v_after)
    integrator.u[2] = v_after

    auto_dt_reset!(integrator)
    set_proposed_dt!(integrator, integrator.dt)
end

cb = ContinuousCallback(condition,nothing,affect_neg! = affect_neg)

u0 = [h0, 0.0]
tspan = (0.0±0.0, 2.0)
prob = ODEProblem(f,u0,tspan)
sol  = solve(prob,Tsit5(),tolerance=1e-6,saveat=0.01,callback=cb)

println("\n... Analytic solution upto second bounce")
t_ev1        = sqrt(2*h0/g)
h_ev1        = -g/2*t_ev1^2 + h0
v_ev1_before = -g*t_ev1
v_ev1_after  = -cof*v_ev1_before

t_ev2        = t_ev1 +2*v_ev1_after/g
h_ev2        = -g/2*(t_ev2-t_ev1)^2 + v_ev1_after*(t_ev2-t_ev1)
v_ev2_before = -v_ev1_after
v_ev2_after  = -cof*v_ev2_before

println("Event at time = ", t_ev1, ", h = ", h_ev1)
println("   v_before = ", v_ev1_before)
println("   v_after  = ", v_ev1_after)

println("Event at time = ", t_ev2, ", h = ", h_ev2)
println("   v_before = ", v_ev2_before)
println("   v_after  = ", v_ev2_after)


# Plott result
tt = Measurements.value.(sol.t)

figure(1)
clf()

h = getindex.(sol.u, 1)
hm = Measurements.value.(h)
hu = Measurements.uncertainty.(h)
hmax = hm + hu
hmin = hm - hu
plot(tt, hmax, label="\$h_{max} in [m]\$")
plot(tt, hmin, label="\$h_{min} in [m]\$")
plot(tt, hm, label="\$h_{mean} in [m]\$")
grid(true)
legend()
title("Bouncing ball with Measurement{Float64}")


figure(2)
clf()

v = getindex.(sol.u, 2)
vm = Measurements.value.(v)
vu = Measurements.uncertainty.(v)
vmax = vm + vu
vmin = vm - vu
plot(tt, vmax, label="\$v_{max} in [m/s]\$")
plot(tt, vmin, label="\$v_{min} in [m/s]\$")
plot(tt, vm, label="\$v_{mean} in [m/s]\$")
grid(true)
legend()
title("Bouncing ball with Measurement{Float64}")

figure(3)
clf()
plot(tt, Measurements.uncertainty.(sol.t), label="uncertainty(time) in [s]")
grid(true)
legend()
title("Bouncing ball with Measurement{Float64}")
end