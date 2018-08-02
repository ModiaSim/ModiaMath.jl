module test_Plot1

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

import ModiaMath

result = Dict{Symbol,AbstractVector}()
t = linspace(0.0, 10.0, 100)
result[:time] = t
result[:phi]  = sin.(t)

ModiaMath.plot(result, :phi, heading="Sine(time)")

end