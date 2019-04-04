# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.Variables (ModiaMath/Variables/_module.jl)
#

# Define ModelVariables struct and functions that holds all information about the variables in a simulation model

isArray(v::ModiaMath.AbstractVariable)     = typeof(v.value) <: AbstractArray
isScalar(v::ModiaMath.AbstractVariable)    = !(typeof(v.value) <: AbstractArray)
isReal(v::ModiaMath.AbstractVariable)      = typeof(v) <: ModiaMath.AbstractRealVariable
valueLength(v::ModiaMath.AbstractVariable) = (isArray(v) ? length(v.value) : 1)
isUsedInAnalysis(v::ModiaMath.AbstractVariable, analysis::AnalysisType) = Int(v.analysis) <= Int(analysis)


#------------------------------- Extract variables from a hierarchical component -------------------------------
function get_ModelVariables_aux!(model, modelType, var, dict, analysis::AnalysisType)
    dict[model] = true   # Mark that model is going to be inspected
    #println("... inspect modelType = ", modelType)

    for i = 1:fieldcount(modelType)
        field = getfield(model, fieldname(modelType, i))
        ftype = typeof(field)
        #println("... fieldname = ", fieldname(modelType,i), ", ftype = ", ftype)
        
        if ftype <: ModiaMath.AbstractVariable 
            if isUsedInAnalysis(field, analysis)
                # Use field since defined for the desired analysis
                push!(var, field)
            end
        elseif ftype <: ModiaMath.AbstractComponentWithVariables && !haskey(dict, field)  # Only inspect further, if model/field was not yet inspected
            get_ModelVariables_aux!(field, ftype, var, dict, analysis)
        end
    end
end


"""
    var = get_ModelVariables(model) - Return all variables defined in model component

The function returns a vector var that contains all AbstractVariables defined in (the hierarchical) model
"""
function get_ModelVariables(model::Any, analysis::AnalysisType)
    var = ModiaMath.AbstractVariable[]
    modelType = typeof(model)

    dict = IdDict{Any,Bool}()
    
    if modelType <: ModiaMath.AbstractComponentWithVariables
        get_ModelVariables_aux!(model, modelType, var, dict, analysis)
    end
    #for v in var
    #   println(v._internal.name)
    #end
    return var   
end


#------------------------------- Construct a variable table for printing/debugging -------------------------------

numericType(v, analysis) = (v.numericType == XD_EXP &&  analysis != ModiaMath.DynamicAnalysis && 
                              Int(v.derivative.analysis) > Int(analysis) ) ? XA : v.numericType 

const staticArrays    = "StaticArrays."
const lenStaticArrays = length(staticArrays)

function shortenedTypeof(value)
    vtype = string(typeof(value))
    if length(vtype) >= lenStaticArrays  && vtype[1:lenStaticArrays ] == staticArrays
        vtype = vtype[lenStaticArrays + 1:end]
    end
    return Symbol(vtype)
end

vecIndex(v) = (isScalar(v) || v.ivar == 0) ? v.ivar : (v.ivar:v.ivar + length(v.value) - 1)

function get_variableTable(var::Vector{ModiaMath.AbstractVariable})
    v_table = DataFrames.DataFrame(name=Symbol[], ValueType=Symbol[], unit=String[],
                numericType=NumericType[], vec=Symbol[], vecIndex=Any[], resultIndex=Any[], 
                der=Symbol[], causality=Causality[], min=Any[], max=Any[], nominal=Any[], fixed=Bool[], 
                start=Any[], info=Symbol[])
   
    for v in var 
        if isReal(v)
            unit        = v.unit
            numericType = v.numericType
            vec         = NumericTypeToVector[Int(v.numericType)]
            der         = typeof(v.derivative) == Nothing ? Symbol("") : instanceName(v.derivative)
            vmin        = v.min
            vmax        = v.max
            vnominal    = v.nominal
        else
            unit        = Symbol("")
            numericType = NoNumericType
            vec         = Symbol("")
            der         = Symbol("")
            vmin        = Symbol("")
            vmax        = Symbol("")
            vnominal    = Symbol("")
        end 
        resultIndex = (isScalar(v) || v.iresult == 0) ? v.iresult : (v.iresult:v.iresult + length(v.value) - 1)
        push!(v_table, [instanceName(v), shortenedTypeof(v.value), unit, numericType, vec, 
                      vecIndex(v), resultIndex, der, v.causality, vmin, vmax, vnominal, v.fixed, v.start, Symbol(v.info)])
    end

    return v_table
end


function indexToString(name, A, linearIndex)
    index = CartesianIndices(A)[linearIndex]
    s     = string(name, "[")
    for i in 1:length(index)
        s = string(s, index[i])
    
        if i == length(index)
            s = string(s, "]")
        else
            s = string(s, ",")
        end
    end
    return s
end


function add_xName!(v::RealVariable, vnumType, xNames, ix_beg)
    if vnumType == XD_EXP || vnumType == XD_IMP || vnumType == XA || vnumType == MUE
        name = string(instanceName(v))
    elseif vnumType == LAMBDA
        name = "integral(" * string(instanceName(v)) * ")"
    else
        error("... should not occur")
    end
    
    if isScalar(v)
        xNames[ix_beg] = string(name)
    else
        for i = 1:length(v.value)
            xNames[ix_beg + i - 1] = indexToString(name, v.value, i)
        end
    end
end


function pushNames!(v::RealVariable, vnumType, xNames, residueNames)
    if vnumType == XD_EXP || vnumType == XD_IMP || vnumType == XA || vnumType == MUE
        push!(xNames, string(instanceName(v)))
    elseif vnumType == LAMBDA
        push!(xNames, "integral(" * string(instanceName(v)) * ")")
    elseif vnumType == FD_IMP || vnumType == FC
        push!(residueNames, string(instanceName(v)))
    end
    
    if vnumType == XD_EXP
        if typeof(v.derivative) == Nothing
            error(instanceName(v), " is defined as XD_EXP, but no variable is defined as its derivative.")
        end
        push!(residueNames, string(instanceName(v.derivative)) * " - derx[.]")
    end
end


"""
    vars = ModiaMath.ModelVariables(model::ModiaMath.AbstractComponentWithVariables; 
                                    analysis::AnalysisType = ModiaMath.DynamicAnalysis)

Return a struct that contains all the variables in `model` in a form so that
fast copying from the integrator interface to the variables can be performed,
and vice versa.
"""
mutable struct ModelVariables
   # Variables used by simulator
    var::Vector{ModiaMath.AbstractVariable}

    # Dimensions of x and derx vector
    nx::Int              # nx = nx_exp + nx_imp + nx_alg + nx_lambda + nx_mue
    nx_exp::Int
    nx_imp::Int
    nx_alg::Int
    nx_lambda::Int
    nx_mue::Int

    # Dimensions of residue vector
    nfd::Int
    nfd_imp::Int    # nfd_exp = nx_exp
    nfc::Int

    # Dimensions of auxiliary variables
    nwr::Int
    nwc::Int

    # Info to copy variables to x/derx/residue/result-values
    x_var::Vector{RealVariable}            # [x_exp_var     , x_imp_var , x_alg_var]
    derx_var::Vector{RealVariable}         # [der(x_imp_var), lambda_var, mue_var]
    residue_var::Vector{RealVariable}      # [fd_imp_var    , fc_var]
    result_var::Vector{ModiaMath.AbstractVariable}   # [time, x_var, derx_var, wr_var, wc_var]
    result_names::Vector{String}           # Names of the (scalarized) result-vector elements
    x_names::Vector{String}                # Names of the (scalarized) x-vector elements

    # Dimensions of the parts of the variable vectors (equivalent to the dimensions of the value vectors like "x" above)
    nx_exp_var::Int
    nx_imp_var::Int
    nx_alg_var::Int
    nx_lambda_var::Int
    nx_mue_var::Int
    nfd_imp_var::Int
    nfc_var::Int
    nwr_var::Int
    nwc_var::Int

    # Dummy
    dummyDifferentialEquation::Bool   # = true, if only dummy differential equation

    function ModelVariables(model::ModiaMath.AbstractComponentWithVariables; analysis::AnalysisType=ModiaMath.DynamicAnalysis)
        # Extract variables from model
        var = get_ModelVariables(model, analysis)
        @assert(length(var) > 0)

        # Determine dimensions of the various parts of the vectors
        ndim = fill(0, 12)
        nvar = fill(0, 12)
        
        for v in var
            numType = Int(numericType(v, analysis))
            ndim[numType] += valueLength(v)
            nvar[numType] += 1
        end

        (nx_exp, nx_imp, nx_alg, nx_lambda, nx_mue, nderx_exp, nderx_imp, nfd_imp, nfc, nwr, nwc, dummy1) = ndim
        (nx_exp_var, nx_imp_var, nx_alg_var, nx_lambda_var, nx_mue_var, nderx_exp_var, nderx_imp_var, nfd_imp_var, nfc_var, nwr_var, nwc_var, dummy2) = nvar
        @assert(nx_exp >= nderx_exp)
        @assert(nx_imp == nderx_imp)
        #println("nx_exp = ", nx_exp, ", nx_imp = ", nx_imp, ", nx_alg = ", nx_alg, ", nx_lambda = ", nx_lambda, ", nx_mue = ", nx_mue)
        nx = nx_exp + nx_imp + nx_alg + nx_lambda + nx_mue

        if nx == 0
            # add dummy differential equation der(x) = -x; x(0)=0
            dummyDifferentialEquation = true
            dummy_x    = RealScalar("_dummy_x"   ; fixed=true, numericType=XD_EXP)
            dummy_derx = RealScalar("_dummy_derx";             numericType=DER_XD_EXP, integral=dummy_x)
            
            pushfirst!(var, dummy_derx)
            pushfirst!(var, dummy_x) 

            nx        = 1
            nx_exp    = 1
            nderx_exp = 1
            nx_exp_var    = 1
            nderx_exp_var = 1
        else
            dummyDifferentialEquation = false
        end
        time = RealScalar("time"; unit="s", causality=Independent, numericType=TIME)
        pushfirst!(var, time) 

        # println("\n... v_table = ", get_variableTable(var))

        # Check dimensions
        if nx_exp + nfd_imp + nfc != nx
            xNames       = String[]
            residueNames = String[]
            for i in eachindex(var)
                v = var[i]
                vnumType = numericType(v, analysis)
                pushNames!(v, vnumType, xNames, residueNames)
            end

            error("The number of x-variables (= ", nx, ") is not identical to the number of equations (= ", nx_exp + nfd_imp + nfc, "):\n",
               "x-variables = ", xNames, "\n",
               "residues = ", residueNames)
        end

        # Allocate variable/index vectors to copy x-, derx-values to variables and variables to residues
        x_var     = Array{ModiaMath.AbstractRealVariable}(undef, nx_exp_var + nx_imp_var + nx_alg_var)   # = [x_exp_var, x_imp_var, x_alg_var]
        ix_exp    = 1
        ix_imp    = nx_exp    + 1
        ix_alg    = ix_imp    + nx_imp
        ix_lambda = ix_alg    + nx_alg
        ix_mue    = ix_lambda + nx_lambda
        ix_exp_var    = 1
        ix_imp_var    = nx_exp_var    + 1
        ix_alg_var    = ix_imp_var    + nx_imp_var
        ix_lambda_var = ix_alg_var    + nx_alg_var
        ix_mue_var    = ix_lambda_var + nx_lambda_var
                   
        derx_var     = Array{ModiaMath.AbstractRealVariable}(undef, nx_imp_var + nx_lambda_var + nx_mue_var)  # = [der(x_imp_var), lambda_var, mue_var]
        iderx_imp    = 1
        iderx_lambda = nx_imp + 1
        iderx_mue    = nx_imp + nx_lambda + 1
        iderx_imp_var    = 1
        iderx_lambda_var = nx_imp_var + 1
        iderx_mue_var    = nx_imp_var + nx_lambda_var + 1
                   
        residue_var  = Array{ModiaMath.AbstractRealVariable}(undef, nfd_imp_var + nfc_var)   # = [fd_imp_var, fc_var]
        ifd_imp = 1
        ifc     = nfd_imp + 1 
        ifd_imp_var = 1
        ifc_var     = nfd_imp_var + 1
                   
        result_var   = Array{ModiaMath.AbstractVariable}(undef, 1 + length(x_var) + nderx_exp_var + length(derx_var) + nwr_var + nwc_var - (dummyDifferentialEquation ? 2 : 0))    # [time, x_var, derx_var, wr_var, wc_var]
        result_names = Array{String}(undef, 1 + nx_exp + nx_imp + nx_alg + nderx_exp + nx_imp + nx_lambda + nx_mue + nwr + nwc - (dummyDifferentialEquation ? 2 : 0))
        result_names[1] = "time"
        iresult     = 1
        iresult_var = 1

        x_names      = Array{String}(undef, nx)

        #println("... length(result_var) = ", length(result_var), ", length(result_names) = ", length(result_names), ", nx_exp = ", nx_exp)

        for i in eachindex(var)
            v      = var[i]
            vnumType = numericType(v, analysis) 
            if vnumType == FD_IMP
                v.ivar                   = nx_exp + ifd_imp
                residue_var[ifd_imp_var] = v
                ifd_imp                 += valueLength(v)
                ifd_imp_var             += 1
            elseif vnumType == FC
                v.ivar               = nx_exp + ifc
                residue_var[ifc_var] = v
                ifc                 += valueLength(v)
                ifc_var             += 1
            else
                if !(dummyDifferentialEquation && (vnumType == XD_EXP || vnumType == DER_XD_EXP))
                    # Determine information about result vectors
                    var_name = string(instanceName(v))
                    if typeof(v.value) == Float64
                        result_names[iresult] = var_name
                        v.iresult = iresult
                        iresult  += 1
                    else
                        v_value = v.value
                        for i in eachindex(v_value)
                            result_names[iresult + i - 1] = indexToString(var_name, v_value, i)
                        end
                        v.iresult = iresult
                        iresult  += length(v_value)
                    end
                    result_var[iresult_var] = v
                    iresult_var += 1
                end

                if vnumType == XD_EXP
                    add_xName!(v, vnumType, x_names, ix_exp)
                    v.ivar             = ix_exp
                    x_var[ix_exp_var]  = v
                    ix_exp            += valueLength(v)
                    ix_exp_var        += 1
                elseif vnumType == XD_IMP
                    add_xName!(v, vnumType, x_names, ix_imp)
                    v.ivar            = ix_imp
                    x_var[ix_imp_var] = v
                    ix_imp           += valueLength(v)
                    ix_imp_var       += 1
                elseif vnumType == XA
                    add_xName!(v, vnumType, x_names, ix_alg)
                    v.ivar            = ix_alg
                    x_var[ix_alg_var] = v
                    ix_alg           += valueLength(v)
                    ix_alg_var       += 1
                elseif vnumType == LAMBDA
                    add_xName!(v, vnumType, x_names, ix_lambda)
                    v.ivar                     = ix_lambda
                    derx_var[iderx_lambda_var] = v
                    ix_lambda        += valueLength(v)
                    iderx_lambda     += valueLength(v)
                    ix_lambda_var    += 1
                    iderx_lambda_var += 1
                elseif vnumType == MUE
                    add_xName!(v, vnumType, x_names, ix_mue)
                    v.ivar                  = ix_mue
                    derx_var[iderx_mue_var] = v
                    ix_mue        += valueLength(v)
                    iderx_mue     += valueLength(v)
                    ix_mue_var    += 1
                    iderx_mue_var += 1
                end
            end
        end

        for v in var
            if v.numericType == DER_XD_IMP
                v.ivar = (v.integral).ivar
                derx_var[iderx_imp_var] = v
                iderx_imp_var += 1   
            elseif v.numericType == DER_XD_EXP
                v.ivar = (v.integral).ivar
            end
        end

        new(var, nx, nx_exp, nx_imp, nx_alg, nx_lambda, nx_mue, nx_exp + nfd_imp, nfd_imp, nfc, nwr, nwc,
          x_var, derx_var, residue_var, result_var, result_names, x_names,
          nx_exp_var, nx_imp_var, nx_alg_var, nx_lambda_var, nx_mue_var, nx_exp_var + nfd_imp_var, nfc_var, nwr_var, nwc_var,
          dummyDifferentialEquation)
    end
end


"""
   table = ModiaMath.get_xTable(vars::ModelVariables)

Function returns a DataFrames tables of all the variables stored in x-vector in `vars`.
"""
function get_xTable(m::ModelVariables)
    x_table = DataFrames.DataFrame(x=Symbol[], name=Symbol[], fixed=Bool[], start=Union{Float64,AbstractArray}[])

    for v in m.x_var
        push!(x_table, [Symbol("x[", vecIndex(v), "]"), instanceName(v), v.fixed, v.start])
    end
    return x_table
end


"""
   table = ModiaMath.get_copyToVariableTable(vars::ModelVariables)

Function returns a DataFrames tables of all the variables in `vars`
that are copied from x/derx to the variables.
"""
function get_copyToVariableTable(m::ModelVariables)
    copyToVariable_table = DataFrames.DataFrame(source=Symbol[], target=Symbol[])

    for v in m.x_var
        push!(copyToVariable_table, [Symbol("x[", vecIndex(v), "]"), instanceName(v)])
    end

    for v in m.derx_var
        push!(copyToVariable_table, [Symbol("derx[", vecIndex(v), "]"), instanceName(v)])
    end
    return copyToVariable_table
end


"""
   table = ModiaMath.get_copyToResidueTable(vars::ModelVariables)

Function returns a DataFrames tables of all the variables in `vars`
that are copied from variables to the residue vector.
"""
function get_copyToResidueTable(m::ModelVariables)
    r_table = DataFrames.DataFrame(source=Symbol[], target=Symbol[])
    x_var = m.x_var
    for i = 1:m.nx_exp_var
        v     = x_var[i]
        der_v = v.derivative
        push!(r_table, [Symbol("derx[", vecIndex(v), "] - ", instanceName(der_v)), Symbol("residue[", vecIndex(v), "]") ])
    end

    for v in m.residue_var
        push!(r_table, [instanceName(v), Symbol("residue[", vecIndex(v), "]") ])
    end
    return r_table
end


"""
   table = ModiaMath.get_copyToResultTable(vars::ModelVariables)

Function returns a DataFrames tables of all the variables in `vars`
that are copied from variables to the result.
"""
function get_copyToResultTable(m::ModelVariables)
    r_table = DataFrames.DataFrame(source=Symbol[], target=Symbol[], start=Union{Float64,AbstractArray}[])

    for v in m.result_var
        resultIndex = isScalar(v) ? v.iresult : (v.iresult:v.iresult + length(v.value) - 1)
        push!(r_table, [instanceName(v), Symbol("result[", resultIndex, "]"), v.start])
    end
    return r_table
end


"""
   table = print_ModelVariables(obj)

Print all the variables in `obj` in form of DataFrames tables.
`obj` can be of type ModiaMath.ModelVariables or ModiaMath.SimulationModel.
"""
function print_ModelVariables(m::ModelVariables)
    variabletable = ModiaMath.get_variableTable(m.var)
    x_table = ModiaMath.get_xTable(m)
    copyToVariableTable = ModiaMath.get_copyToVariableTable(m)
    copyToResidueTable = ModiaMath.get_copyToResidueTable(m)
    copyToResultTable = ModiaMath.get_copyToResultTable(m)

    print("\n\nvariables: ");                show(variabletable      , splitcols=true, summary=false)
    print("\n\n\nx vector: ");               show(x_table            , splitcols=true, summary=false)
    print("\n\n\ncopy to variables: ");      show(copyToVariableTable, splitcols=true, summary=false)
    print("\n\n\ncopy to residue vector: "); show(copyToResidueTable , splitcols=true, summary=false)
    print("\n\n\ncopy to results: ");        show(copyToResultTable  , splitcols=true, summary=false)
    print("\n")
end


"""
    copy_x_and_derx_to_variables!(time::Float64, x::Vector{Float64}, 
                                  derx::Vector{Float64}, vars::ModelVariables)

Copy `x`and `derx`of the integrator interface to the model variables `vars`.
"""
function copy_x_and_derx_to_variables!(time::Float64, x::Vector{Float64}, derx::Vector{Float64}, m::ModelVariables)
    @assert(length(x)    == m.nx)
    @assert(length(derx) == m.nx)

    m.var[1].value = time
    for v in m.x_var
        #println("... typeof(", fullName(v), ") = ", typeof(v), ", isimmutable(v) = ", isimmutable(v))
        if isScalar(v)
            v.value = x[ v.ivar ]
        elseif isimmutable(v.value)
            # v is an immutable array (e.g. SVector{3,Float64})
            v.value = x[ v.ivar:v.ivar + length(v.value) - 1 ]
        else
            vv = v.value
            for j in 1:length(vv)
                vv[j] = x[ v.ivar + j - 1 ]
            end
        end
    end

    for v in m.derx_var
        if isScalar(v)
            v.value = derx[ v.ivar ]
        elseif isimmutable(v.value)
            # v is an immutable array (e.g. SVector{3,Float64})
            v.value = derx[ v.ivar:v.ivar + length(v.value) - 1 ]
        else
            vv = v.value
            for j in 1:length(vv)
                vv[j] = derx[ v.ivar + j - 1 ]
            end
        end
    end

    return nothing
end


"""
    copy_start_to_x!(vars::ModelVariables, x::Vector{Float64}, x_fixed::Vector{Bool}, x_nominal::Vector{Float64})

Copy `start`, `fixed` and `nominal`values of variables `vars` to `x`, `x_fixed`, and `x_nominal` vectors.
"""
function copy_start_to_x!(m::ModelVariables, x::Vector{Float64}, x_fixed::Vector{Bool}, x_nominal::Vector{Float64})
    @assert(length(x)         == m.nx)
    @assert(length(x_fixed)   == m.nx)
    @assert(length(x_nominal) == m.nx) 
  
    for v in m.x_var
        ibeg = v.ivar
        if isScalar(v)
            x[ibeg]         = v.start
            x_fixed[ibeg]   = v.fixed
            x_nominal[ibeg] = v.nominal
        else
            vv = v.start
            for j in 1:length(vv)
                k = ibeg + j - 1
                x[ k ]       = vv[j]
                x_fixed[   k ] = v.fixed 
                x_nominal[ k ] = v.nominal 
            end
        end
    end
    return nothing
end
copy_start_to_x!(m::ModelVariables, x::Vector{Float64}, x_fixed::Vector{Bool}) = copy_start_to_x!(m, x, x_fixed, fill(1.0,m.nx))


"""
    copy_variables_to_residue!(vars::ModelVariables, x::Vector{Float64},
                               derx::Vector{Float64}, residue::Vector{Float64})

Copy variables `vars` to `residue` vector and include the inherent residue equations
of `XD_EXP` variables (`residue = der(v) - derx`).
"""
function copy_variables_to_residue!(m::ModelVariables, x::Vector{Float64}, derx::Vector{Float64}, residue::Vector{Float64})
    @assert(length(x)       == m.nx)
    @assert(length(derx)    == m.nx)
    @assert(length(residue) == m.nx)

    # If dummy equations present, compute derivative of dummy equation
    if m.dummyDifferentialEquation
        # der(x) = -x
        m.var[3].value = -m.var[2].value
    end

    # Generate residues of explicitly solvable derivatives
    x_var = m.x_var
    for i = 1:m.nx_exp_var
        v = x_var[i]
        if isScalar(v)
            residue[v.ivar] = derx[v.ivar] - v.derivative.value
        else
            der_vv = v.derivative.value
            for j in 1:length(der_vv)
                k = v.ivar + j - 1
                residue[ k ] = derx[ k ] - der_vv[j]
            end            
        end
    end

    # Copy residue variables of fd_imp and fc to residue vector
    for v in m.residue_var
        if isScalar(v)
            residue[ v.ivar ] = v.value
        else
            vv = v.value
            for j in 1:length(vv)
                residue[ v.ivar + j - 1 ] = vv[j]
            end
        end
    end

    return nothing
end



function get_variableValueTable(m::ModelVariables)
    v_table = DataFrames.DataFrame(simulator=Symbol[], variable=Symbol[], value=Any[])

    for v in m.x_var
        push!(v_table, [Symbol("x[", vecIndex(v), "]"), instanceName(v), v.value])
    end

    for v in m.derx_var
        push!(v_table, [Symbol("derx[", vecIndex(v), "]"), instanceName(v), v.value])
    end
    return v_table
end


function get_residueValueTable(m::ModelVariables, derx_integrator::Vector{Float64})
    r_table = DataFrames.DataFrame(simulator=Symbol[], variable=Symbol[], value=Any[])
    x_var = m.x_var
    for i = 1:m.nx_exp_var
        v     = x_var[i]
        der_v = v.derivative
        value = derx_integrator[ vecIndex(v) ] - der_v.value
        push!(r_table, [Symbol("residue[", vecIndex(v), "]"), Symbol("derx[", vecIndex(v), "] - ", instanceName(der_v)), value ])
    end

    for v in m.residue_var
        push!(r_table, [Symbol("residue[", vecIndex(v), "]"), instanceName(v), v.value ])
    end
    return r_table
end


