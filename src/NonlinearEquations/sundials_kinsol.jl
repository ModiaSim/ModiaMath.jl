# License for this file: MIT (expat)
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

import Sundials
import ModiaMath
using LinearAlgebra

@noinline function old_cfunction(f, r, a)
    ccall(:jl_function_ptr, Ptr{Cvoid}, (Any, Any, Any), f, r, a)
end


mutable struct NonlinearEquationsInfo
    extraInfo                # Model-specific extra information
    name::String             # Name of equation system
    ny::Int                  # Number of equations (length of y-vector)
    getResidues!::Function   # Function of the nonlinear equation system
    y0::Vector{Float64}      # The initial y vector (to print in error messages)
    lastNorm_r::Float64
    lastrScaledNorm_r::Float64
    kin_mem::Ptr{Cvoid}      # KINSOL pointer (to access all KINgetXXX functions)

    function NonlinearEquationsInfo(name::String, ny::Int, getResidues!::Function; extraInfo=nothing)
        @assert(ny >= 0)
        new(extraInfo, name, ny, getResidues!, zeros(ny), 1.0, 1.0)
    end
end

function kinsol_f(y::Sundials.N_Vector, r::Sundials.N_Vector, eqInfo::NonlinearEquationsInfo)
    eqInfo.getResidues!(eqInfo, Sundials.asarray(y), Sundials.asarray(r))
    return Cint(0)   # indicates normal return
end

@static if VERSION < v"0.7.0-DEV.2005"
    const kinsol_fc = cfunction(kinsol_f, Cint, (Sundials.N_Vector, Sundials.N_Vector, Ref{NonlinearEquationsInfo}))
end

function kinsol_ErrHandlerFn(error_code::Cint, KINmodule::Cstring, KINfunction::Cstring,
                             message::Cstring, eqInfo::NonlinearEquationsInfo)
    if error_code > 0
        # Print information text
        println("\n\n!!! Warning from ModiaMath.NonlinearEquations.KINSOL: ", unsafe_string(KINfunction),
               "(...) returned with a [", unsafe_string(KINmodule), "] error_code  = ",
               error_code, ":\n", unsafe_string(message), "\n")
    else
        simulationModel = eqInfo.extraInfo
        tables = ModiaMath.getVariableAndResidueValues(simulationModel)

        if tables != nothing
            println("\n\nLast used values in model:\n\n",
                  tables[1], "\n\n",
                  tables[2], "\n\n")
        end

        if error_code == -11
            str1 = "\nIt might be that the Jacobian is singular (= there are redundant equations).\n"
        else
            str1 = "\n"
        end

        if typeof(simulationModel) <: ModiaMath.AbstractSimulationModel
            simState = simulationModel.simulationState
            x_names     = String[simState.getVariableName(simState.model, ModiaMath.DAE.Category_X, i)  for i = 1:simState.nx]

            error("\n\n!!! Error from ModiaMath.NonlinearEquations.KINSOL: ", unsafe_string(KINfunction),
                  "(...) returned with a [", unsafe_string(KINmodule), "] error:\n    ",
                  unsafe_string(message), "\nModiaMath info:\nlastNorm(r) = ", eqInfo.lastNorm_r,
                  ", lastNorm(rScaled*r) = ", eqInfo.lastrScaledNorm_r, ".", str1,
                  string(simState.name), ": time = " , string(simState.time) ,
                 ", stepsize of implicit Euler step = " , string(simState.hev) ,
                 ", scaleConstraintsAtEvents = " , string(simState.scaleConstraintsAtEvents) ,
                 "\nx_names     = " , string(x_names) ,
                 "\nx_start     = " , string(eqInfo.y0) ,
                 "\nx_fixed     = " , string(simState.x_fixed) ,
                 "\nx_nominal   = " , string(simState.x_nominal) ,
                 "\nx_yScale    = " , string(simState.yScale) ,
                 "\nx_rScale    = " , string(simState.rScale) ,
                 "\nx           = " , string(simState.xev) ,
                 "\nderx        = " , string(simState.derxev) ,
                 "\nresidues    = " , string(simState.residues) ,
                 "\nnx          = " , string(simState.nx) ,
                 "\nnd          = " , string(simState.nd) ,
                 "\nnc          = " , string(simState.nc) ,
                 "\nnw          = " , string(simState.nw) ,
                 "\nnz          = " , string(simState.nz) ,
                 "\ntolerance   = " , string(simState.tolerance) ,
                 "\nFTOL        = " , string(simState.FTOL))

        else
            error("\n\n!!! Error from ModiaMath.NonlinearEquations.KINSOL: ", unsafe_string(KINfunction),
                  "(...) returned with a [", unsafe_string(KINmodule), "] error:\n    ",
                  unsafe_string(message), "\nModiaMath info:\nlastNorm(r) = ", eqInfo.lastNorm_r,
                  ", lastNorm(rScaled*r) = ", eqInfo.lastrScaledNorm_r, ".", str1)
        end
    end
    return nothing
end

@static if VERSION < v"0.7.0-DEV.2005"
    const kinsol_ErrHandlerFnc = cfunction(kinsol_ErrHandlerFn, Void, (Cint, Cstring, Cstring, Cstring, Ref{NonlinearEquationsInfo}))
end


function solveNonlinearEquations!(eqInfo::NonlinearEquationsInfo, y::Vector{Float64};
                                  FTOL::Float64=eps(Float64)^(1 / 3),
                                  yScale::Vector{Float64}=ones(length(y)),
                                  rScale::Vector{Float64}=ones(length(y)))
    # Create KINSOL and info structure
    @assert(length(y) == eqInfo.ny)
    kmem = Sundials.KINCreate()
    if kmem == C_NULL
        error("ModiaMath.NonlinearEquations.KINSOL.solveNonlinearEquations!: Failed to allocate KINSOL solver object")
    end
    eqInfo.kin_mem = kmem
    eqInfo.y0 .= y

    # Run KINSOL
    try
        # Set error handler function
        @static if VERSION >= v"0.7.0-DEV.2005"
            Sundials.KINSetErrHandlerFn(kmem, old_cfunction(kinsol_ErrHandlerFn, Nothing, Tuple{Cint,Cstring,Cstring,Cstring,Ref{typeof(NonlinearEquationsInfo)}}), pointer_from_objref(eqInfo))
        else
            Sundials.KINSetErrHandlerFn(kmem, kinsol_ErrHandlerFnc, pointer_from_objref(eqInfo))
        end

        # Initialize KINSOL
        @static if VERSION >= v"0.7.0-DEV.2005"
            Sundials.KINInit(kmem, old_cfunction(kinsol_f, Cint, Tuple{Sundials.N_Vector,Sundials.N_Vector,Ref{typeof(NonlinearEquationsInfo)}}), Sundials.NVector(y))
        else
            Sundials.KINInit(kmem, kinsol_fc, y)
        end
        Sundials.KINSetUserData(kmem, eqInfo)
        Sundials.KINSetFuncNormTol(kmem, FTOL)
        Sundials.KINSetScaledStepTol(kmem, FTOL * FTOL)

        # Set maximum allowable scaled length mxnewtstep of the Newton step
        # KINSOL defines a default of mxnewtstep = 1000*norm(y.*yScale).
        # This fails if y is a zero vector. Therefore, mxnewtstep is changed to
        # take only yScale into account.
        Sundials.KINSetMaxNewtonStep(kmem, 1000.0 * norm(yScale, Inf))

        # Set linear solver
        #   nnz_jac = nnz(sim.jac)
        #   IDAKLU(mem, ny, nnz_jac)
        #   IDAKLUSetOrdering(mem, KLUorderingChoice)
        #   IDASlsSetSparseJacFn(mem);
        #else
        A = Sundials.SUNDenseMatrix(length(y), length(y))
        LS = Sundials.SUNDenseLinearSolver(Sundials.NVector(y), A)
        Sundials.KINDlsSetLinearSolver(kmem, LS, A)
        #end
        strategy = Sundials.KIN_LINESEARCH

        @static if VERSION >= v"0.7.0-DEV.2005"
            Sundials.KINSol(kmem, Sundials.NVector(y), strategy, Sundials.NVector(yScale), Sundials.NVector(rScale))
        else
            Sundials.KINSol(kmem, y, strategy, yScale, rScale)
        end

    finally
        # Free allocated memory
        Sundials.KINFree([kmem])
    end

    return nothing
end

end