# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module Simulate_PendulumDAE

Simulate PendulumDAE model (= index-3 model transformed to index-1 model)
"""
module Simulate_PendulumDAE

include(joinpath("models", "PendulumDAE.jl"))
import .PendulumDAE

import ModiaMath

model  = PendulumDAE.Model()
result = ModiaMath.simulate!(model, stopTime=2.0, log=true) 

ModiaMath.plot(result, [(:x, :y), (:vx, :vy), :lambda, :mue], heading="Simulate_PendulumDAE.jl (index3 reduced to index1)")

end