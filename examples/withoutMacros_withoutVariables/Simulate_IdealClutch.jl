# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module Simulate_IdealClutch

Simulate models/IdealClutch.jl model.
    
Demonstrates/tests varying number of constraint equations fc leading to Dirac impulses at events
"""
module Simulate_IdealClutch

include(joinpath("models", "IdealClutch.jl"))
import .IdealClutch
import ModiaMath

model  = IdealClutch.Model()
result = ModiaMath.simulate!(model, stopTime=500.0; log=true) 
ModiaMath.plot(result, ("w1", "w2", "C.v"), heading="Simulate_IdealClutch.jl")

end