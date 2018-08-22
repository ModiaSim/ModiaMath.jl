# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module Simulate_PendulumODE2

Simulate PendulumODE model.
"""
module Simulate_PendulumODE2

include(joinpath("models", "PendulumODE2.jl"))
import .PendulumODE2
import ModiaMath

model  = PendulumODE2.Model(L=0.8, m=0.5, d=0.0)
result = ModiaMath.simulate!(model, stopTime=5.0) 

heading = "Simulate_PendulumODE2.jl (ODE as index-0 DAE)"
ModiaMath.plot(result, [:phi, :w, :der_w], figure=1, heading=heading)
ModiaMath.plot(result, :ry, xAxis=:rx, figure=2, heading=heading)

end