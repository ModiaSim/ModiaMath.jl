# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module ModiaMath.NonlinearEquations

This module contains functions to solve nonlinear algebraic equations:

- [`solveOneNonlinearEquation`](@ref):
  Function that computes the solution of one non-linear algebraic equation `y=f(u)`
  in one unknown `u` in a reliable and efficient way using Brents algorithm.

- [`KINSOL`](@ref):
  Module containing functions to solve a system of nonlinear systems of equations
  with Sundials KINSOL. The module is designed so that the same system is
  solved several times with KINSOL as needed by Modia simulations
  (auxiliary memory is only allocated once and not for every call).
  KINSOL is used in ModiaMath to solve nonlinear algebraic equations during
  initialization and at events of a simulation.

# Main developer
[Martin Otter](https://rmc.dlr.de/sr/de/staff/martin.otter/),
[DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)
"""
module NonlinearEquations

export KINSOL, NonlinearEquationsInfo, solveNonlinearEquations!
export solveOneNonlinearEquation

# using/imports
import ModiaMath

# include code
include("solveOneNonlinearEquation.jl")
include("sundials_kinsol.jl")

end
