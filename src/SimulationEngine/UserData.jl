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


function sol_f!(simModel::IntegratorData, sim, t::Float64, _y::Vector{Float64}, _yp::Vector{Float64}, r::Vector{Float64}, hcur)
#    simModel.hcur[1] = simModel.sol_mem.hcur
#    simModel.order[1] = simModel.sol_mem.order
    stat           = simModel.statistics
    stat.hMin      = min(stat.hMin, simModel.hcur[1])
    stat.hMax      = max(stat.hMax, simModel.hcur[1])
    stat.orderMax  = max(stat.orderMax, simModel.order[1])
    stat.nResidues += 1

    #simModel.y = _y

    #simModel.yp = _yp

    sim                  = simModel.simulationState
    sim.time             = t
    simModel.last_t      = t
    simModel.last_norm_y = norm(simModel.y, Inf)

    # Check that simModel.y is still identical to _y


    ModiaMath.DAE.getResidues!(simModel.model, sim, t, _y, _yp, r, hcur)
    #println(_y, _yp, t, r, hcur )
    return nothing
end


#------- root finding
function sol_g(t::Float64, y::Vector{Float64}, simModel::IntegratorData)
    sim      = simModel.simulationState
    sim.time = t
    yp = simModel.yp
    ModiaMath.DAE.getEventIndicators!(simModel.model, sim, t, y, yp, simModel.z)
    return simModel.z   # indicates normal return
end


#------- full jacobian
function idasol_fulljac(t::Float64, cj::Float64, _y::Vector{Float64}, _yp::Vector{Float64}, _r::Vector{Float64},
                        _fulljac::Array{Float64}, simModel::IntegratorData,
                        tmp1::Vector{Float64}, tmp2::Vector{Float64}, tmp3::Vector{Float64})
    # Compute full Jacobian
    sim      = simModel.simulationState
    sim.time = t
    simModel.hcur[1] = simModel.sol_data.hcur
    simModel.eweight = simModel.sol_data.weights


    simModel.y = _y

    simModel.yp = _yp

    r = _r
    ModiaMath.DAE.computeJacobian!(simModel.model, sim, t, y, yp, r, simModel.fulljac, simModel.hcur[1], cj, simModel.eweight)
    simModel.statistics += sim.nx

    return 0   # indicates normal return
end




function find_roots(t::Array{Float64, 1}, ys::Array{Array{Float64, 1}, 1}, yps::Array{Array{Float64,1},1}, simModel::IntegratorData,  resprob!, abstol::Float64, reltol::Array{Float64,1})
    #solve again with better precision to be sure
    FTOL=eps(Float64)^(1 / 3)
    sim = simModel.simulationState
    nsol_f!(r, y) = ModiaMath.DAE.getResidues!(simModel.model, sim, t[1], ys[1], yps[1], r, simModel.hcur[1])
    nsol = nlsolve(nsol_f!, ys[1], ftol=FTOL)
    y = nsol.zero
    prob = DAEProblem{true}(resprob!, yps[1], y, (t[1], t[end]), differential_vars=[true,true,false])
    sol = solve(prob, dassl(), reltol = reltol, abstol = abstol/100)


    return sol.t[end], sol.u[end], sol.du[end]
end
