module test_Plot4

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
result["time"] = t * u"s"
result["phi"]  = sin.(t)u"rad"
result["phi2"] = 0.5 * sin.(t)u"rad"
result["w"]    = cos.(t)u"rad/s"
result["w2"]   = 0.6 * cos.(t)u"rad/s"
result["r"]    = hcat(0.4 * cos.(t), 0.5 * sin.(t), 0.3*cos.(t))
result["r2"]   = hcat(0.1*cos.(t), 0.2*cos.(t), 0.3*cos.(t), 0.4*cos.(t), 0.5*cos.(t),
                      0.6*cos.(t), 0.7*cos.(t), 0.8*cos.(t), 0.9*cos.(t), 1.0*cos.(t), 1.1*cos.(t))

ModiaMath.plot(result, [ (:phi), (:phi, :phi2, :w), (:w, :w2) ], heading="Vector of plots")

ModiaMath.plot(result, [ "r", "r2[3]", "r2"], heading="Vector of plots", figure=2)

end