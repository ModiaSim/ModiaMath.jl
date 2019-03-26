# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module
#   ModiaMath.NonlinearEquations (ModiaMath/NonlinearEquations/_module.jl)
#

"""
    module KINSOL - Solve nonlinear equation system with Sundials KINSOL

The goal is to solve the same system several times with KINSOL.
"""

module KINSOL

import Flux
import ModiaMath
using LinearAlgebra
using Flux

mutable struct NonlinearEquationsInfo
    extraInfo                # Model-specific extra information
    name::String             # Name of equation system
    ny::Int                  # Number of equations (length of y-vector)
    getResidues!::Function   # Function of the nonlinear equation system
    y0::Vector{Float64}      # The initial y vector (to print in error messages)
    lastNorm_r::Float64
    lastrScaledNorm_r::Float64

    function NonlinearEquationsInfo(name::String, ny::Int, getResidues!::Function; extraInfo=nothing)
        @assert(ny >= 0)
        new(extraInfo, name, ny, getResidues!, zeros(ny), 1.0, 1.0)
    end
end


function solveNonlinearEquations!(eqInfo::NonlinearEquationsInfo, y::Vector{Float64})
                                  #inputs are the same for backwards compatibility

    function j!(J, x)
        J = Flux.Tracker.jacobian(eqInfo.getResidues!, x)
    end
        # Display the ODE with the initial parameter values.
    print("nonlinear")
    nlsolve(eqInfo.getResidues!, j!, y)

    return nothing
end

end
