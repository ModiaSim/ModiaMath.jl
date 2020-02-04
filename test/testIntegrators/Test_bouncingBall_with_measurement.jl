module Test_bouncingBall_with_measurement

#=
Strange result of bouncing ball with Measurement.jl

https://github.com/JuliaPhysics/Measurements.jl/issues/64

Issue fixed in Measurements version 2.2.0
=#

using DifferentialEquations, Measurements, PyPlot


const g   = 9.81       # Gravity constant
const h0  = 1.0 ± 0.1  # Initial height
const cor = 0.7 ± 0.05  # Coefficient of restitution


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
    v_after = -cor*v_before
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
v_ev1_after  = -cor*v_ev1_before
v_ev1_0_46   = -g*(0.46-t_ev1) + v_ev1_after

t_ev2        = t_ev1 +2*v_ev1_after/g
h_ev2        = -g/2*(t_ev2-t_ev1)^2 + v_ev1_after*(t_ev2-t_ev1)
v_ev2_before = -v_ev1_after
v_ev2_after  = -cor*v_ev2_before

println("Event at time = ", t_ev1, ", h = ", h_ev1)
println("   v_before = ", v_ev1_before)
println("   v_after  = ", v_ev1_after)
println("   v(0.46)  = ", v_ev1_0_46)

println("Event at time = ", t_ev2, ", h = ", h_ev2)
println("   v_before = ", v_ev2_before)
println("   v_after  = ", v_ev2_after)


# Plot result
tm = Measurements.value.(sol.t)
tu = Measurements.uncertainty.(sol.t)

function plotUncertaintyOfEvents(tm, tu, vm, vu; color=[0.0, 0.0, 1.0, 0.5], linewidth=3.0, color_uncertainty=[0, 0.16, 0.78, 0.2])
    first = true
    for i in 1:length(tm)-1
        if tm[i] == tm[i+1] && tu[i] > 0.0
            tt = [tm[i]-tu[i], tm[i]+tu[i]]
            vlow = vm[i] - vu[i]
            vup  = vm[i+1] + vu[i+1]
            fill_between(tt,[vlow, vlow], [vup, vup], color=color_uncertainty)
            if first
                first = false
                plot(tt, [0.0, 0.0], color=color, linewidth=linewidth, label="uncertainty(time)")
            else
                plot(tt, [0.0, 0.0], color=color, linewidth=linewidth)
            end
        end
    end
end

figure(1)
clf()

h = getindex.(sol.u, 1)
hm = Measurements.value.(h)
hu = Measurements.uncertainty.(h)
hmax = hm + hu
hmin = hm - hu
#plot(tt, hmax, label="\$h_{max} in [m]\$")
#plot(tt, hmin, label="\$h_{min} in [m]\$")

# Remove uncertainty values with h < 0
for i in 1:length(hmin)
    if hmin[i] < 0
        hmin[i] = 0
    end
end
fill_between(tm,hmin,hmax, color=[0, 0.16, 0.78, 0.2], label="uncertainty(h)")
plot(tm, hm, color=[0.0, 0.0, 1.0, 1.0], label="\$h_{mean}\$ in [m]")
grid(true)
plotUncertaintyOfEvents(tm,tu,hm,hu)
xlabel("time in [s]")
legend()
text(0.25, 1.1, "\$h_{0} = \$" * string(h0))
text(0.25, 1.05, "\$cor = \$" * string(cor))
title("Bouncing ball with Measurement{Float64}")


figure(2)
clf()

v = getindex.(sol.u, 2)
vm = Measurements.value.(v)
vu = Measurements.uncertainty.(v)
vmax = vm + vu
vmin = vm - vu
#plot(tm, vmax, label="\$v_{max}\$")
#plot(tm, vmin, label="\$v_{min}\$")
fill_between(tm,vmin,vmax, color=[0, 0.16, 0.78, 0.2], label="uncertainty(v)")
plot(tm, vm, color=[0.0, 0.0, 1.0, 1.0], label="\$v_{mean}\$ in [m/s]")
grid(true)
plotUncertaintyOfEvents(tm,tu,vm,vu)
xlabel("time in [s]")
legend(loc="lower right")
text(0, 3.5, "\$h_{0} = \$" * string(h0))
text(0, 3.0, "\$cor = \$" * string(cor))
title("Bouncing ball with Measurement{Float64}")

figure(3)
clf()
plot(tm, tu, label="uncertainty(time) in [s]")
grid(true)
legend()
title("Bouncing ball with Measurement{Float64}")
end