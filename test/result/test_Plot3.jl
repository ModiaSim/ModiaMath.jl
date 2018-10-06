module test_Plot3

import ModiaMath

#  Desired:
#    using Test
#    using Unitful
#  
#  In order that these packages need not to be defined in the user environment, they are included via ModiaMath:

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
    t = linspace(0.0, 10.0, 100)
else
    using ModiaMath.Test
    t = range(0.0, stop=10.0, length=100)
end

using ModiaMath.Unitful

result = Dict{AbstractString,Any}()
result["time"]  = t
result["phi"]   = sin.(t)u"rad"
result["phi2"]  = 0.5 * sin.(t)u"rad"
result["w"]     = cos.(t)u"rad/s"

println("\n... Next plot should give a warning:")
ModiaMath.plot(result, (:phi, :phi2, :w, :signalNotDefined), heading="Sine(time)", figure=1)

ModiaMath.plot(result, [:phi, :phi2, :w], heading="Sine(time)", figure=2)
ModiaMath.plot(result, :phi, xAxis=:w, heading="phi=f(w)", figure=3)

println("\n... Next plot should give a warning:")
ModiaMath.plot(result, :phi, xAxis=:xAxisNotDefined, heading="phi=f(w)", figure=4)

# Add new simulation result
result["phi"]  = 1.2*result["phi"]
result["phi2"] = 1.1*result["phi2"]
result["w"]    = 0.5*result["w"]
ModiaMath.plot(result, (:phi, :phi2, :w), figure=1, prefix="Sim 2: ", reuse=true)
ModiaMath.plot(result, [:phi, :phi2, :w], figure=2, prefix="Sim 2: ", reuse=true)

end