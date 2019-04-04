# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module Simulate_PendulumDAE2

Simulate PendulumDAE model with x_fixed=true.
"""
module Simulate_PendulumDAE2

include(joinpath("models", "PendulumDAE.jl"))
import .PendulumDAE

import ModiaMath

model  = PendulumDAE.Model(x_fixed=true)
result = ModiaMath.simulate!(model, stopTime=2.0, log=true) 

ModiaMath.plot(result, [(:x, :y), (:vx, :vy), :lambda, :mue], heading="Simulate_PendulumDAE2.jl (index3 reduced to index1)")

end