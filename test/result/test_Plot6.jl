module test_Plot6

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


# Define signals
result = Dict{AbstractString,Any}("time" => t      , "phi1" => sin.(t)     , "phi2" => 0.5 * sin.(t),
                                  "w1"   => cos.(t), "w2"   => 0.6 * cos.(t))

# Plots
ModiaMath.plot(result, :phi1)                                        # 1 signal in one diagram
ModiaMath.plot(result, (:phi1, :phi2, :w1), figure=2)                # 3 signals in one diagram
ModiaMath.plot(result, [:phi1, :phi2, :w1], figure=3)                # 3 diagrams in form of a vector (every diagram has one signal) 
ModiaMath.plot(result, ["phi1" "phi2";
                        "w1"   "w2"   ], figure=4)                   # 4 diagrams in form of a matrix (every diagram has one signal)
ModiaMath.plot(result, [ (:phi1, :phi2), (:w1) ], figure=5)           # 2 diagrams in form of a vector
ModiaMath.plot(result, [ (:phi1,)           (:phi2, :w1);
                         (:phi1, :phi2, :w1)  (:w2,)     ], figure=6)   # 4 diagrams in form of a matrix
ModiaMath.plot(result, :w1, xAxis=:phi1, figure=7)                   # Plot w1=f(phi1) in one diagram

end