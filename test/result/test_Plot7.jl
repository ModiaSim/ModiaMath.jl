module test_Plot7

import ModiaMath

# Desired:
#   using Test
#
# In order that Test needs not to be defined in the user environment, it is included via ModiaMath:
using ModiaMath.Test
t = range(0.0, stop=10.0, length=100)



# Define signals
mutable struct MyThermodynamicState
    p::Float64
    T::Float64
end

specificEnthalpy(state::MyThermodynamicState) = 2.0*state.T
dynamicViscosity(state::MyThermodynamicState) = 2.0*state.p

const dependentVariables = Dict{AbstractString,Function}("h"   => specificEnthalpy,
                                                         "eta" => dynamicViscosity)

ModiaMath.variablesDependingOnStruct(state::MyThermodynamicState) = dependentVariables


state = MyThermodynamicState[]
for i = 1:length(t)
    push!(state, MyThermodynamicState(i*2.0,i*3.0))
end

result = Dict{AbstractString,Any}("time" => t, "state" => state)

# Plots
ModiaMath.plot(result, ("state.p", "state.T", "state.h", "state.eta"), figure=7)

end