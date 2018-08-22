# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    modul Simulate_BouncingBall

Simulate bouncing ball (model with state events and discontinuous change of DAE states).
"""
module Simulate_BouncingBall

include(joinpath("models", "BouncingBall.jl"))
import .BouncingBall
import ModiaMath

# Simulate
model  = BouncingBall.Model()
result = ModiaMath.simulate!(model, stopTime=3.0, log=true) 

ModiaMath.plot(result, [:h, (:v, :flying)], heading="Simulate_BouncingBall.jl")

end