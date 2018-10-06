module test_Plot2

import ModiaMath

# Desired:
#   using Test
#
# In order that Test needs not to be defined in the user environment, it is included via ModiaMath:
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
    t = linspace(0.0, 10.0, 100)
else
    using ModiaMath.Test
    t = range(0.0, stop=10.0, length=100)
end


series = Dict{AbstractString,Any}()
series["time"] = t
series["phi"]  = sin.(t)
series["w"]    = cos.(t)
series["r"]    = hcat(0.4 * cos.(t), 0.5 * sin.(t), 0.3*cos.(t))

time = ModiaMath.RealScalar(:time, unit="s")
phi  = ModiaMath.RealScalar(:phi, unit="rad")
w    = ModiaMath.RealScalar(:w, unit="rad/s")
r    = ModiaMath.RealSVector3(:r, unit="m")
var  = Dict{Symbol,Any}(:time => time, :phi => phi, :w => w, :r => r)

result = ModiaMath.ResultWithVariables(series, var, "ModiaMath/test/Result/test_Plot2.jl")

ModiaMath.plot(result, :phi    , prefix="sim 1: ", heading="Sine(time)")
ModiaMath.plot(result, "r"     , prefix="sim 1: ", figure=2)
ModiaMath.plot(result, "r[2]"  , figure=3)

println("\n... Next plot should give a warning:")
ModiaMath.plot(result, "r[2:3]", figure=4)

println("\n... Next plot should give a warning:")
ModiaMath.plot(result, "phi", xAxis="r", figure=4)


# Add next simulation run to plot
result.series["phi"] = 1.2*sin.(t)
result.series["r"]   = hcat(0.2 * cos.(t), 0.3 * sin.(t), 0.5*cos.(t))
ModiaMath.plot(result, :phi, heading="Sine(time)", figure=1, prefix="sim 2: ", reuse=true)
ModiaMath.plot(result, :r  ,                       figure=2, prefix="sim 2: ", reuse=true)


end