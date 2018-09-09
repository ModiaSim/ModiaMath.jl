module test_Plot3

import ModiaMath
using  ModiaMath.Unitful


@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
    t = linspace(0.0, 10.0, 100)
else
    using ModiaMath.Test
    t = range(0.0, stop=10.0, length=100)
end
 



@testset "ModiaMath.Result: test Result" begin

    result = Dict{Symbol,AbstractVector}()
    result[:time]  = t
    result[:phi]   = sin.(t)u"rad"
    result[:phi2]  = 0.5 * sin.(t)u"rad"
    result[:w]     = cos.(t)u"rad/s"

    ModiaMath.plot(result, (:phi, :phi2, :w, :signalNotDefined), heading="Sine(time)", figure=1)
    ModiaMath.plot(result, [:phi, :phi2, :w], heading="Sine(time)", figure=2)
    ModiaMath.plot(result, :phi, xAxis=:w, heading="phi=f(w)", figure=3)
    ModiaMath.plot(result, :phi, xAxis=:xAxisNotDefined, heading="phi=f(w)", figure=4)
end

end