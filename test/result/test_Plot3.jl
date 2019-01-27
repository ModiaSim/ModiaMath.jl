module test_Plot3

import ModiaMath

#  Desired:
#    using Test
#    using Unitful
#  
#  In order that these packages need not to be defined in the user environment, they are included via ModiaMath:

using ModiaMath.Test
t = range(0.0, stop=10.0, length=100)


using ModiaMath.Unitful

result = Dict{AbstractString,Any}()
result["time"]    = t
result["phi"]     = sin.(t)u"rad"
result["phi2"]    = 0.5 * sin.(t)u"rad"
result["w"]       = cos.(t)u"rad/s"
result["phi_max"]     = 1.1u"rad"
result["phi_max_int"] = 1
result["open"]        = false

println("\n... Next plot should give a warning:")
ModiaMath.plot(result, (:phi_max, :phi_max_int, :open, :phi, :phi2, :w, :signalNotDefined), heading="Sine(time)", figure=1)

ModiaMath.plot(result, [:phi, :phi2, :w], heading="Sine(time)", figure=2)
ModiaMath.plot(result, :phi, xAxis=:w, heading="phi=f(w)", figure=3)

println("\n... Next plot should give a warning:")
ModiaMath.plot(result, :phi, xAxis=:xAxisNotDefined, heading="phi=f(w)", figure=4)

# Print result variables
println("\n... result variables = ", ModiaMath.resultTable(result))



# Add new simulation result
result["phi"]  = 1.2*result["phi"]
result["phi2"] = 1.1*result["phi2"]
result["w"]    = 0.5*result["w"]
ModiaMath.plot(result, (:phi, :phi2, :w), figure=1, prefix="Sim 2: ", reuse=true)
ModiaMath.plot(result, [:phi, :phi2, :w], figure=2, prefix="Sim 2: ", reuse=true)

end