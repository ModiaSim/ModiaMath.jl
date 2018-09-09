module test_Plot4

import ModiaMath
using  ModiaMath.Unitful


@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
    t = linspace(0.0, 10.0, 100)
else
    using ModiaMath.Test
    t = range(0.0, stop=10.0, length=100)
end



result = Dict{Symbol,AbstractVector}()
result[:time] = t * u"s"
result[:phi]  = sin.(t)u"rad"
result[:phi2] = 0.5 * sin.(t)u"rad"
result[:w]    = cos.(t)u"rad/s"
result[:w2]   = 0.6 * cos.(t)u"rad/s"

ModiaMath.plot(result, [ (:phi), (:phi, :phi2, :w), (:w, :w2) ], heading="Vector of plots")

end