################################################
#
# Simulate model\PendulumODE
#
# Author: Martin Otter, DLR-SR 
#         (first version: Nov. 19, 2016)
#
################################################

"""
    Simulate_PendulumODE2 - Simulate PendulumODE model
"""
module Simulate_PendulumODE2

include(joinpath("models","PendulumODE2.jl"))
import .PendulumODE2
import ModiaMath

model  = PendulumODE2.Model(L=0.8, m=0.5, d=0.0)
result = ModiaMath.simulate!(model, stopTime=5.0) 

heading="simulationWithoutMacro/Simulate_PendulumODE2.jl"
ModiaMath.plot(result, [:phi, :w, :der_w], figure=1, heading=heading)
ModiaMath.plot(result, :ry, xAxis=:rx    , figure=2, heading=heading)

end