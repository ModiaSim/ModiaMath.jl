# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module ModiaMath.NonlinearEquations

Solve nonlinear algebraic equations

# Main developer
Martin Otter, [DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)
"""
module NonlinearEquations

export KINSOL, NonlinearEquationsInfo, solveNonlinearEquations!

# using/imports
import ModiaMath

# include code
include("sundials_kinsol.jl")

end
