# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module
#   ModiaMath.SimulationEngine(ModiaMath/SimulationEngine/_module.jl)
#

# Methods that involve IDA UserData


#=
"KLU sparse matrix type"
type SlsMat
   M::Cint                        # number of rows
   N::Cint                        # number of columns
   NNZ::Cint                      # maximum number of nonzero entries in the matrix
   data::Ptr{Sundials.realtype}   # data[NNZ]    - value of nonzero entries
   rowvals::Ptr{Cint}             # rowvals[NNZ] - row indices of data (index starts at zero)
   colptrs::Ptr{Cint}             # colptrs[N+1] - index of the first column entry into the
                                  #                data and rowvals arrays (index starts at zero)
                                  #                The last entry contains the total number of nonzero values
                                  #                in the matrix and hence points one past the end
                                  #                of the active data in the data and rowvals arrays.
end

const SlsMat_Ptr = Ptr{SlsMat}


"""
    CSCtoSlsMat!(mat, matKLU)

Copy SparseMatrixCSC mat to SlsMat data structure matKLU
(storage for both objects are provided in the calling program)
"""
function CSCtoSlsMat!(mat::SparseMatrixCSC{Float64,Cint}, matKLU::SlsMat)
   for i = 1:length(mat.nzval)
      unsafe_store!(matKLU.data   , mat.nzval[i], i)
      unsafe_store!(matKLU.rowvals, mat.rowval[i]-1, i)
   end

   for i = 1:length(mat.colptr)
      unsafe_store!(matKLU.colptrs, mat.colptr[i]-1, i)
   end
end
=#



##################################################################
#
# Wrappers to IDA code
#
##################################################################


function idasol_f(t::Sundials.realtype, _y::Sundials.N_Vector, _yp::Sundials.N_Vector, _r::Sundials.N_Vector, simModel::IntegratorData)
    Sundials.IDAGetCurrentStep(simModel.ida_mem, simModel.hcur)
    Sundials.IDAGetCurrentOrder(simModel.ida_mem, simModel.order)
    stat           = simModel.statistics
    stat.hMin      = min(stat.hMin, simModel.hcur[1])
    stat.hMax      = max(stat.hMax, simModel.hcur[1])
    stat.orderMax  = max(stat.orderMax, simModel.order[1])
    stat.nResidues += 1

    sim                  = simModel.simulationState
    sim.time             = t
    simModel.last_t      = t
    simModel.last_norm_y = norm(simModel.y, Inf)

    # Check that simModel.y is still identical to _y
    if pointer(simModel.y) === Sundials.__N_VGetArrayPointer_Serial(_y)
       # If this assumption is true, no unnecessary memory is allocated via Sundials.asarray(..)
        y  = simModel.y
    else
        y  = Sundials.asarray(_y)
        println("\n!!! Info message from ModiaMath.simulate:\n",
                 "    idasol_f assumption SimModel.y === _y is not valid; using asarray(_y).")
    end

    if pointer(simModel.yp) === Sundials.__N_VGetArrayPointer_Serial(_yp)
       # If this assumption is true, no unnecessary memory is allocated via Sundials.asarray(..)
        yp = simModel.yp
    else
        yp = Sundials.asarray(_yp)
        println("\n!!! Info message from ModiaMath.simulate:\n",
                 "    idasol_f assumption SimModel.yp === _yp is not valid; using asarray(_yp).")
    end

    ModiaMath.DAE.getResidues!(simModel.model, sim, t, y, yp, simModel.r, simModel.hcur[1])

    # Copy simModel.r to _r
    unsafe_copyto!(Sundials.__N_VGetArrayPointer_Serial(_r), pointer(simModel.r), simModel.simulationState.nx)
    return Cint(0)   # indicates normal return
end


@noinline function old_cfunction(f, r, a)
    ccall(:jl_function_ptr, Ptr{Cvoid}, (Any, Any, Any), f, r, a)
end


#------- root finding
function idasol_g(t::Sundials.realtype, y::Sundials.N_Vector, yp::Sundials.N_Vector, gout::Ptr{Sundials.realtype}, simModel::IntegratorData)
    simModel.statistics.nZeroCrossings += 1
    sim      = simModel.simulationState
    sim.time = t
    ModiaMath.DAE.getEventIndicators!(simModel.model, sim, t, Sundials.asarray(y), Sundials.asarray(yp), simModel.z)

    for i = 1:simModel.simulationState.nz
        unsafe_store!(gout, simModel.z[i], i)
    end
    return Cint(0)   # indicates normal return
end


#------- full jacobian
function idasol_fulljac(t::Sundials.realtype, cj::Sundials.realtype, _y::Sundials.N_Vector, _yp::Sundials.N_Vector, _r::Sundials.N_Vector,
                        _fulljac::Sundials.SUNMatrix, simModel::IntegratorData,
                        tmp1::Sundials.N_Vector, tmp2::Sundials.N_Vector, tmp3::Sundials.N_Vector)
    # Compute full Jacobian
    sim      = simModel.simulationState
    sim.time = t
    Sundials.IDAGetCurrentStep(simModel.ida_mem, simModel.hcur)
    IDAGetErrWeights(simModel.ida_mem, simModel.eweight)


    # Check that simModel.y is still identical to _y
    if pointer(simModel.y) === Sundials.__N_VGetArrayPointer_Serial(_y)
       # If this assumption is true, no unnecessary memory is allocated via Sundials.asarray(..)
        y  = simModel.y
    else
        y  = Sundials.asarray(_y)
        println("\n!!! Info message from ModiaMath.simulate:\n",
                 "    idasol_fulljac assumption SimModel.y === _y is not valid; using asarray(_y).")
    end

    if pointer(simModel.yp) === Sundials.__N_VGetArrayPointer_Serial(_yp)
       # If this assumption is true, no unnecessary memory is allocated via Sundials.asarray(..)
        yp = simModel.yp
    else
        yp = Sundials.asarray(_yp)
        println("\n!!! Info message from ModiaMath.simulate:\n",
                 "    idasol_fulljac assumption SimModel.yp === _yp is not valid; using asarray(_yp).")
    end

    r = Sundials.asarray(_r)
    ModiaMath.DAE.computeJacobian!(simModel.model, sim, t, y, yp, r, simModel.fulljac, simModel.hcur[1], cj, simModel.eweight)
    simModel.statistics += sim.nx

    # Copy simModel.fulljac to fulljac
    unsafe_copyto!(Sundials.__N_VGetArrayPointer_Serial(_fulljac), pointer(simModel.fulljac), sim.nx)
    return Cint(0)   # indicates normal return
end




#=
#-------- Jacobian
function idasol_sjac(t::Sundials.realtype, c_j::Sundials.realtype, y::Sundials.N_Vector, yp::Sundials.N_Vector, r::Sundials.N_Vector,
                     jacKLU::SlsMat, simModel::IntegratorData,
                     tmp1::Sundials.N_Vector, tmp2::Sundials.N_Vector, tmp3::Sundials.N_Vector)
    # Compute non-zero values of jac
    sim      = simModel.simulationState
    sim.time = t
    h        = IDAGetCurrentStep(simModel.ida_mem)
    IDAGetErrWeights(simModel.ida_mem, simModel.eweight)

    ModiaMath.SparseJacobian.computeJacobian!(ModiaMath.DAE.getResidues!, simModel.model, simModel.simulationState,
                                              t, Sundials.asarray(y), Sundials.asarray(yp), Sundials.asarray(r),
                                              sim.cg, c_j, h, Sundials.asarray(simModel.eweight), Sundials.asarray(tmp1), sim.jac)
    simModel.statistics.nResidues += sim.cg.ngroups

    # Copy userData.jac to jacKLU
    #   function call gives a Julia crash: CSCtoSlsMat!(simModel.cm.jac, jacKLU)
    #   Manually inlining the function works
    jac = sim.jac
    for i = 1:length(jac.nzval)
      unsafe_store!(jacKLU.data   , jac.nzval[i], i)
      unsafe_store!(jacKLU.rowvals, jac.rowval[i]-1, i)
    end

    for i = 1:length(jac.colptr)
      unsafe_store!(jacKLU.colptrs, jac.colptr[i]-1, i)
    end

    return Cint(0)   # indicates normal return
end


const idasol_sjacc = cfunction(idasol_sjac, Int32, (Sundials.realtype, Sundials.realtype, Sundials.N_Vector, Sundials.N_Vector, Sundials.N_Vector, Ref{SlsMat},
                                                    Ref{IntegratorData}, Sundials.N_Vector, Sundials.N_Vector, Sundials.N_Vector))
=#


function idasol_ErrHandlerFn(error_code::Cint, IDAmodule::Cstring, IDAfunction::Cstring,
                             message::Cstring, simModel::IntegratorData)
    modelName = simModel.simulationState.name
    if error_code > 0
        # Print information text
        println("\n\n!!! Warning from ModiaMath.simulate(", modelName, ", ...): ", unsafe_string(IDAfunction),
               "(...) returned with [", unsafe_string(IDAmodule), "] error_code  = ",
               error_code, ":\n", unsafe_string(message), "\n")
    elseif error_code < 0
        # Raise error
        time   = simModel.last_t
        norm_y = simModel.last_norm_y

        if time > typemin(Float64)
            error("\n\n!!! Error from ModiaMath.simulate(", modelName, ", ...): ", unsafe_string(IDAfunction),
                "(...) returned with an [", unsafe_string(IDAmodule), "] error:\n",
                unsafe_string(message), "\n( norm( x(t=", @sprintf("%0.3g",time), " s) ) = ", @sprintf("%0.3g",norm_y),
                "; if this value is large, the model is unstable;\n",
                "  if message 'mxstep steps taken before reachint tout', set maxNumberOfSteps to a value > 500 )\n")
        else
            error("\n\n!!! Error from ModiaMath.simulate(", modelName, ", ...): ", unsafe_string(IDAfunction),
                "(...) returned with an [", unsafe_string(IDAmodule), "] error:\n",
                unsafe_string(message), "\n")
        end
    end
    nothing
end
