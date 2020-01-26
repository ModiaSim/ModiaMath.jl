module Test_pendulum_dae

#=
SavingCallback for DAEs: Not clear how to get du #547

https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/547

=#


using DifferentialEquations, Measurements, StaticArrays, PyPlot, Sundials

mutable struct Model{R<:Real}
    # parameters
    m::R  # mass
    L::R  # length
    d::R  # damping
    g::R  # gravity

	# variables
    ϕ::R  # angle
    ω::R  # der(angle)
    α::R  # der(der(angle))
	residue::R
    r::SVector{2,R}   # absolute position of end point

    der_x::Vector{R}

    function Model{R}(m::R, L::R, d::R, g::R) where {R<:Real}
        new(m, L, d, g, 0, 0, 0, 0, SVector{2,R}(0,0))
    end
end

Model(;m::R=0, L::R=0, d::R=0, g::R=0) where {R<:Real} = Model{R}(m,L,d,g)


"simplePendulum: Compute m::Model at actual time instant"
function pendulum!(m::Model, t, x, der_x)::Nothing
    m.ϕ = x[1]
    m.ω = x[2]
    m.α = (-m.m*m.g*m.L*sin(m.ϕ) - m.d*m.ω) / (m.m*m.L^2)
	m.residue = der_x[2] - m.α
    m.r = @SVector [m.L*sin(m.ϕ), -m.L*cos(m.ϕ)]
    return nothing
end

"SimplePendulum1: Called by integrator"
function pendulum1!(residue, der_x, x, m::Model, t)
    pendulum!(m, t, x, der_x)
    residue[1] = der_x[1] - m.ω
    residue[2] = m.residue
end

"SimplePendulum2: Called at every communication point"
function pendulum2!(x, t, integrator)
    m::Model = integrator.p
    # tmps = get_tmp_cache(integrator)
    # println("typeof(tmps) = ", typeof(tmps), ", length(tmps) = ", length(tmps))
    if t==0
        m.der_x = copy(integrator.du)   # Since IDA gives an error for integrator(t, Val{1]}) at the initial time instant
    else
        integrator(m.der_x, t, Val{1})
    end
    pendulum!(m, t, x, m.der_x)
    return deepcopy(m)
end

t_start = 0.0
t_inc   = 0.01
t_end   = 6.0
tspan   = (t_start, t_end)
tspan2  = t_start:t_inc:t_end


# Problem 1 with Float64
model1 = Model(m=1.0, L=1.0, d=0.2, g=9.81)
println("model1 = $model1")
x1₀     = [π/2, 0.0]

    # Initialize
    der_x1₀ = [0.0, 0.0]
    pendulum!(model1, 0.0, x1₀, der_x1₀)
    der_x1₀[1] = model1.ω
    der_x1₀[2] = model1.α

saved_values1=SavedValues(typeof(x1₀[1]), Model)
cb1 = SavingCallback(pendulum2!, saved_values1, saveat=tspan2)

prob1 = DAEProblem(pendulum1!, der_x1₀, x1₀, tspan, model1)
sol1 = solve(prob1, IDA(), reltol = 1e-6, saveat=t_inc, callback=cb1, save_timeseries=false)

#@time sol1 = solve(prob1, IDA(), reltol = 1e-6, saveat=t_inc)
#@time sol1 = solve(prob1, IDA(), reltol = 1e-6, saveat=t_inc)


# switch off storing of x, der_x in sol1 (but keep in saved_values1)
# sol1 = solve(prob1, IDA(), reltol = 1e-6, saveat=t_end, callback=cb1, save_timeseries=false)

t = saved_values1.t
m = saved_values1.saveval

figure(1)
clf()
plot(sol1.t, getindex.(sol1.u, 1), label="\$\\varphi\$")
#ϕ = getfield.(m, :ϕ)
#plot(t, ϕ, label="\$\\varphi\$")
grid(true)
legend()
title("Pendulum with Float64")

figure(2)
clf()

r = getfield.(m, :r)
r1 = getindex.(r,1)
r2 = getindex.(r,2)
plot(t, r1, label="r1")
plot(t, r2, label="r2")
grid(true)
legend()
title("Pendulum with Float64")

figure(3)
clf()
residue = getfield.(m, :residue)
plot(t, residue, label="residue")
grid(true)
legend()
title("Pendulum with Float64")


end