module test_Plot3

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

import ModiaMath
using Unitful

result = Dict{Symbol,AbstractVector}()
t = linspace(0.0, 10.0, 100)
result[:time] = t
result[:phi]   = sin.(t)u"rad"
result[:phi2]  = 0.5*sin.(t)u"rad"
result[:w]     = cos.(t)u"rad/s"

ModiaMath.plot(result, (:phi, :phi2, :w), heading="Sine(time)", figure=1)
ModiaMath.plot(result, [:phi, :phi2, :w], heading="Sine(time)", figure=2)
ModiaMath.plot(result, :phi, xAxis=:w, heading="phi=f(w)", figure=3)

end