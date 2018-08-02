################################################
#
# Author: Martin Otter, DLR-SR 
#         (first version: Nov. 26, 2016)
#
################################################

"""
    Simulate_BouncingBall - Simulate bouncing ball
"""
module Simulate_BouncingBall

include(joinpath("models","BouncingBall.jl"))
import .BouncingBall
import ModiaMath

# Simulate
model  = BouncingBall.Model()
result = ModiaMath.simulate!(model, stopTime=3.0, log=true) 

ModiaMath.plot(result, [:h, (:v, :flying)], heading="Simulate_BouncingBall.jl")

end