module Test_ode_pendulum_with_Interval

#=
Same as Test_ode_pendulum_with_Measurement.jl but using IntervalArithmetic instead of Measurement

Result: Very long simulation time + result is useless (too large interval)
=#

using DifferentialEquations, IntervalArithmetic, StaticArrays, PyPlot

Base.Float64(x::IntervalArithmetic.Interval{Float64}) = (x.hi + x.lo)/2

mutable struct Model{R<:Real, F<:Real}
    m::R  # mass
    L::R  # length
    d::R  # damping
    g::F  # gravity (not used as Measurment)
    ϕ::R  # angle
    ω::R  # der(angle)
    α::R  # der(der(angle))
    r::SVector{2,R}   # absolute position of end point

    function Model{R,F}(m::R, L::R, d::R, g::F) where {R<:Real, F<:Real}
        new(m, L, d, g, 0, 0, 0, SVector{2,R}(0,0))
    end
end

Model(;m::R=0, L::R=0, d::R=0, g::F=0) where {R<:Real, F<:Real} = Model{R,F}(m,L,d,g)

#Model(m::R, L::R, d::R, g::F) where {R<:Real, F<:Real} = Model{R,F}(m,L,d,g)
#Model(;m::Real=0, L::Real=0, d::Real=0, g::Real=0) = Model(promote(m,L,d),...)


"simplePendulum: Compute m::Model at actual time instant"
function simplePendulum!(m::Model, t, x)::Nothing
    m.ϕ = x[1]
    m.ω = x[2]
    m.α = (-m.m*m.g*m.L*sin(m.ϕ) - m.d*m.ω) / (m.m*m.L^2)
    m.r = @SVector [m.L*sin(m.ϕ), -m.L*cos(m.ϕ)]
    return nothing
end

"SimplePendulum1: Called by integrator"
function simplePendulum1!(der_x, x, m::Model, t)
    simplePendulum!(m, t, x)
    der_x[1] = m.ω
    der_x[2] = m.α
end

"SimplePendulum2: Called at every communication point"
function simplePendulum2(x, t, integrator)
    m::Model = integrator.p
    simplePendulum!(m, t, x)
    return deepcopy(m)
end

t_start = 0.0
t_inc   = 0.01
t_end   = 7.0
tspan   = (t_start, t_end)
tspan2  = t_start:t_inc:t_end


# Problem 1 with Float64
model1 = Model(m=1.0, L=1.0, d=0.2, g=9.81)
println("model1 = $model1")
x1₀ = [π/2, 0.0]

saved_values1=SavedValues(typeof(x1₀[1]), Model)
cb1 = SavingCallback(simplePendulum2, saved_values1, saveat=tspan2)

prob1 = ODEProblem(simplePendulum1!, x1₀, tspan, model1)
sol1 = solve(prob1, Tsit5(), reltol = 1e-6, saveat=t_inc, callback=cb1)
prob1b = deepcopy(prob1)
sol3 = solve(prob1, Tsit5(), reltol = 1e-6, saveat=t_inc)

figure(1)
clf()
plot(sol1.t, getindex.(sol1.u, 1), label="\$\\varphi\$")
grid(true)
legend()
title("SimplePendulum with Float64")

figure(2)
clf()
t = saved_values1.t
m = saved_values1.saveval
r = getfield.(m, :r)
r1 = getindex.(r,1)
r2 = getindex.(r,2)
plot(t, r1, label="r1")
plot(t, r2, label="r2")
grid(true)
legend()
title("SimplePendulum with Float64")


# Problem 2 with Measurement
model2 = Model(m=1.0 ± 0.1, L=1.00 ± 0.01, d=0.2± 0.05, g=9.81)
println("model2 = $model2")
println("typeof(model2.r) = ", typeof(model2.r))
x2₀ = [π/2 ± 0, 0.0 ± 0]

saved_values2=SavedValues(typeof(x2₀[1]), Model)
cb2 = SavingCallback(simplePendulum2, saved_values2, saveat=tspan2)

prob2 = ODEProblem(simplePendulum1!, x2₀, tspan, model2)
sol2  = solve(prob2, Tsit5(), reltol = 1e-6, saveat=t_inc, callback=cb2)

ϕ = getindex.(sol2.u, 1)
ϕmax = getfield.(ϕ,:hi)
ϕmin = getfield.(ϕ,:lo)

figure(3)
clf()
plot(sol2.t, ϕmax, label="\$\\varphi_{max}\$")
plot(sol2.t, ϕmin, label="\$\\varphi_{min}\$" )
grid(true)
legend()
title("SimplePendulum with Interval{Float64}")

end