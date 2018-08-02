################################################
#
# Simulate model\PendulumODE
#
# Author: Martin Otter, DLR-SR 
#         (first version: Nov. 19, 2016)
#
################################################

"""
    Simulate_PendulumODE - Simulate PendulumODE model
"""
module Simulate_PendulumODE


include(joinpath("models","PendulumODE.jl"))
import .PendulumODE
import ModiaMath

model  = PendulumODE.Model(L=0.8, m=0.5, d=0.2)
result = ModiaMath.simulate!(model, stopTime=5.0, log=true)

heading="simulationWithoutMacro/Simulate_PendulumODE.jl"
ModiaMath.plot(result, [:phi, :w, :der_w], figure=1, heading=heading)
ModiaMath.plot(result, :ry, xAxis=:rx    , figure=2, heading=heading)

end