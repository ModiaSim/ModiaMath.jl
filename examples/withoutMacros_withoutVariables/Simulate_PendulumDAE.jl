################################################
#
# Simulate model\PendulumODE
#
# Author: Martin Otter, DLR-SR 
#         (first version: Dec. 4, 2016)
#
################################################

"""
    Simulate_PendulumDAE - Simulate PendulumDAE model
"""
module Simulate_PendulumDAE

include(joinpath("models","PendulumDAE.jl"))
import .PendulumDAE

import ModiaMath

model  = PendulumDAE.Model()
result = ModiaMath.simulate!(model, stopTime=2.0, log=true) 

ModiaMath.plot(result, [(:x, :y), (:vx, :vy), :lambda], heading="simulationWithoutMacro/Simulate_PendulumDAE.jl")

end