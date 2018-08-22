# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module Simulate_SimpleStateEvents

Simulate SimpleStateEvents model.
    
Demonstrates/tests simple state events
"""
module Simulate_SimpleStateEvents

include(joinpath("models", "SimpleStateEvents.jl"))
import .SimpleStateEvents
import ModiaMath

model  = SimpleStateEvents.Model()
result = ModiaMath.simulate!(model, stopTime=10.0; log=true) 
ModiaMath.plot(result, (:s, :v, :f), heading="sSimulate_SimpleStateEvents.jl")

end