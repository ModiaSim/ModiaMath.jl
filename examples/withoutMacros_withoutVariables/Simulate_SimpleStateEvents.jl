################################################
#
# Author: Martin Otter, DLR-SR 
#         (first version: Nov. 26, 2016)
#
################################################

"""
    Simulate_SimpleStateEvents - Simulate SimpleStateEvents model
    
Demonstrates/tests simple state events
"""
module Simulate_SimpleStateEvents

include(joinpath("models","SimpleStateEvents.jl"))
import .SimpleStateEvents
import ModiaMath


model  = SimpleStateEvents.Model()
result = ModiaMath.simulate!(model, stopTime=10.0; log=true) 
ModiaMath.plot(result, (:s, :v, :f), heading="simulationWithoutMacro/Simulate_SimpleStateEvents.jl")


end