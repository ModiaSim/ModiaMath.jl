module test_Plot1

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
    t = linspace(0.0, 10.0, 100)
else
    using Test
    t = range(0.0, stop=10.0, length=100)
end

import ModiaMath

result = Dict{Symbol,AbstractVector}()

result[:time] = t
result[:phi]  = sin.(t)

ModiaMath.plot(result, :phi, heading="Sine(time)")

end