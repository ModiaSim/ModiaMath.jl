module test_Plot2

import ModiaMath
using Test

t = range(0.0, stop=10.0, length=100)

series = Dict{Symbol,Any}()
series[:time] = t
series[:phi]  = sin.(t)
series[:w]    = cos.(t)
series[:r]    = hcat(0.4 * cos.(t), 0.5 * sin.(t))

time = ModiaMath.RealScalar(:time, unit="s")
phi  = ModiaMath.RealScalar(:phi, unit="rad")
w    = ModiaMath.RealScalar(:w, unit="rad/s")
r    = ModiaMath.RealSVector3(:r, unit="m")
var  = Dict{Symbol,Any}(:time => time, :phi => phi, :w => w, :r => r)

result = ModiaMath.ResultWithVariables(series, var, "ModiaMath/test/Result/test_Plot2.jl")

ModiaMath.plot(result, :phi, heading="Sine(time)")
ModiaMath.plot(result, :r, figure=2)
ModiaMath.plot(result, "r[2]", figure=3)


end
