# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module Simulate_PendulumODE 

Simulate PendulumODE model.
"""
module Simulate_PendulumODE

include(joinpath("models", "PendulumODE.jl"))
import .PendulumODE
import ModiaMath


model  = PendulumODE.Model(L=0.8, m=0.5, d=0.2)
result = ModiaMath.simulate!(model, stopTime=5.0, log=true)

heading = "Simulate_PendulumODE.jl (ODE as index-0 DAE)"
ModiaMath.plot(result, [:phi, :w, :der_w], figure=1, heading=heading)
ModiaMath.plot(result, :ry, xAxis=:rx, figure=2, heading=heading)

end