module Test_ode_save_func_callback

#=
StackOverFlowError with two usages of save_func #548

https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/548

Issue is fixed by setting:

   ArrayInterface.ismutable(::Type{Model}) = true

=#



using DifferentialEquations, StaticArrays

mutable struct Model{R<:Real}
    m::R  # mass
    L::R  # length
    d::R  # damping
    g::R  # gravity
    ϕ::R  # angle
    ω::R  # der(angle)
    α::R  # der(der(angle))
    r::SVector{2,R}   # absolute position of end point

    function Model{R}(m::R, L::R, d::R, g::R) where {R<:Real}
        new(m, L, d, g, 0, 0, 0, SVector{2,R}(0,0))
    end
end

# Due to StackOverFlowError https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/548
using ArrayInterface
ArrayInterface.ismutable(::Type{Model}) = true


Model(;m::R=0, L::R=0, d::R=0, g::R =0) where {R<:Real} = Model{R}(m,L,d,g)


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


model1 = Model(m=1.0, L=1.0, d=0.2, g=9.81)
x1₀ = [π/2, 0.0]

saved_values1=SavedValues(typeof(x1₀[1]), Model)
cb1 = SavingCallback(simplePendulum2, saved_values1, saveat=tspan2)

prob1 = ODEProblem(simplePendulum1!, x1₀, tspan, model1)
sol1 = solve(prob1, Tsit5(), reltol = 1e-6, saveat=t_inc, callback=cb1)

# Triggers stack-overflow error if ArrayInterface.ismutable(::Type{Model}) = true is not set
sol2 = solve(prob1, Tsit5(), reltol = 1e-6, saveat=t_inc, callback=cb1)

end