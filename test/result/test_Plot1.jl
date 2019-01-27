module test_Plot1

import ModiaMath

# Desired:
#   using Test
#
# In order that Test needs not to be defined in the user environment, it is included via ModiaMath:
using ModiaMath.Test
t = range(0.0, stop=10.0, length=100)




result = Dict{AbstractString,Any}()

result["time"] = t
result["phi"]  = sin.(t)

ModiaMath.plot(result, :phi, heading="Sine(time)")

# Test also that warn works for non-existing variables
println("... Next plot should give a warning:")
ModiaMath.plot(result, :SignalNotDefined, heading="Sine(time)", figure=2)

# Print result variables
println("\n... result variables = ", ModiaMath.resultTable(result))

end