################################################
#
# Simulate model\PendulumDAE with x_fixed=true
#
# Author: Martin Otter, DLR-SR 
#         (first version: Aug. 15, 2017)
#
################################################

"""
    Simulate_PendulumDAE2 - Simulate PendulumDAE model with x_fixed=true
"""
module Simulate_PendulumDAE2

include(joinpath("models","PendulumDAE.jl"))
import .PendulumDAE

import ModiaMath

model  = PendulumDAE.Model(x_fixed=true)
result = ModiaMath.simulate!(model, stopTime=2.0, log=true) 

ModiaMath.plot(result, [(:x, :y), (:vx, :vy), :lambda], heading="simulationWithoutMacro/Simulate_PendulumDAE2.jl")

end