# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
     module Simulate_FreeBodyRotation

Simulate body with unconstrained rotational motion in 3D (using quaternions; = index-1 model).
"""
module Simulate_FreeBodyRotation

include(joinpath("models", "FreeBodyRotation.jl"))
import .FreeBodyRotation
import ModiaMath

model  = FreeBodyRotation.Model()
result = ModiaMath.simulate!(model, stopTime=5.0, log=true) 

ModiaMath.plot(result, [("Q[1]", "Q[2]", "Q[3]", "Q[4]"), 
                        ("w[1]", "w[2]", "w[3]"), 
                        ("der(w[1])", "der(w[2])", "der(w[3])")], heading="simulationWithoutMacro/Simulate_FreeBodyRotation.jl")

end
