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

import NLsolve, LineSearches
import ModiaMath
using LinearAlgebra
using NLsolve


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


function solveNonlinearEquations!(eqInfo::NonlinearEquationsInfo, y::Vector{Float64};
                                  FTOL::Float64=eps(Float64)^(1 / 2),
                                  yScale::Vector{Float64}=ones(length(y)),
                                  rScale::Vector{Float64}=ones(length(y)))
                                  #inputs are the same for backwards compatibility
    eqInfo.y0 .= y
    itnum = Int32(1000.0 * norm(yScale, Inf))
    nsol_f!(r, y) = eqInfo.getResidues!(eqInfo, y, r)
    nsol = nlsolve(nsol_f!, y, method = :newton, ftol=FTOL, iterations=itnum, linesearch=LineSearches.HagerZhang())
    #nsol = nlsolve(nsol_f!, y, ftol=FTOL, iterations=itnum)
    eqInfo.y0 .= nsol.zero
    return nothing
end

end
