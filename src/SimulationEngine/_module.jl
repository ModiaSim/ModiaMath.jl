# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control


"""
    module ModiaMath.SimulationEngine

Simulation engine for implicit index 1 DAE models with events.

# Main developer

[Martin Otter](https://rmc.dlr.de/sr/de/staff/martin.otter/),
[DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)
"""
module SimulationEngine

export simulate!

import Sundials
import ModiaMath
import ModiaMath.DAE
using  StaticArrays
@eval using Printf
using LinearAlgebra


# include code
# include("IDA_UserData.jl")  # is included within simulate.jl
include("simulate.jl")

end
