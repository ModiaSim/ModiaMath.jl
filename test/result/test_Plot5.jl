module test_Plot5

import ModiaMath

# Desired:
#   using Test
#
# In order that Test needs not to be defined in the user environment, it is included via ModiaMath:
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
    t = linspace(0.0, 10.0, 100)
else
    using ModiaMath.Test
    t = range(0.0, stop=10.0, length=100)
end



series = Dict{AbstractString,Any}()
series["time"] = t
series["phi"]  = sin.(t)
series["phi2"] = 0.5 * sin.(t)
series["w"]    = cos.(t)
series["w2"]   = 0.6 * cos.(t)
series["r"]    = hcat(0.5 * sin.(t), 0.2 * sin.(t), 0.3 * cos.(t))


time = ModiaMath.RealScalar(:time, unit="s")
phi  = ModiaMath.RealScalar(:phi,  unit="rad")
phi2 = ModiaMath.RealScalar(:phi2, unit="rad")
w    = ModiaMath.RealScalar(:w,    unit="rad/s")
w2   = ModiaMath.RealScalar(:w2,   unit="rad/s")
r    = ModiaMath.RealSVector3(:r,  unit="m")
var  = Dict{Symbol,Any}(:time => time, :phi => phi, :phi2 => phi2, :w => w, :w2 => w2, :r => r)

result = ModiaMath.ResultWithVariables(series, var, "")

ModiaMath.plot(result, [ (:phi, :r)      (:phi, :phi2, :w);
                         (:w, :w2, :phi2) (:phi, :w)      ], heading="Matrix of plots")

end