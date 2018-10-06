# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.Variables (ModiaMath/Variables/_module.jl)
#

# Define SimulationModel struct and functions


"""
    sm = ModiaMath.SimulationModel(model::ModiaMath.AbstractComponentWithVariables;
                                   startTime = 0.0, stopTime = 1.0, tolerance = 1e-4,
                                   interval = (stopTime-startTime)/500.0)

Return a simulationModel `sm` struct that can be simulated with [`ModiaMath.simulate!`](@ref).
The given `model` is typically constructed with the [`@component`](@ref) macro.
As keyword arguments, default values for `startTime`, `stopTime`, `tolerance`, `interval` can be given.
"""
mutable struct SimulationModel <: ModiaMath.AbstractSimulationModel
    modelName::String
    simulationState::ModiaMath.SimulationState
    var::ModelVariables
    model::ModiaMath.AbstractComponentWithVariables 

    function SimulationModel(model::ModiaMath.AbstractComponentWithVariables;
                             startTime = 0.0,
                             stopTime  = 1.0,
                             tolerance = 1e-4,
                             interval  = (stopTime-startTime)/500.0,
                             hev = 1e-8,
                             scaleConstraintsAtEvents::Bool = true)
        modelName = ModiaMath.componentName(model)
        var = ModelVariables(model)
        x = zeros(var.nx)
        x_fixed = fill(false,var.nx)
        ModiaMath.copy_start_to_x!(var,x,x_fixed)

        simulationState = ModiaMath.SimulationState(modelName, getModelResidues!, x, getVariableName; 
                                                    x_fixed = x_fixed, nc = var.nfc,
                                                    getResultNames = getResultNames, storeResult! = storeVariables!,
                                                    getResult        = getResult,
                                                    defaultStartTime = startTime,
                                                    defaultStopTime  = stopTime,
                                                    defaultTolerance = tolerance,
                                                    defaultInterval  = interval,
                                                    hev = hev,
                                                    scaleConstraintsAtEvents = scaleConstraintsAtEvents)

        # model._internal.simulationState = simulationState
        @static if VERSION >= v"0.7.0-DEV.2005"
            new(String(modelName), simulationState, var, model)
        else
            new(modelName, simulationState, var, model)
        end
    end
end 

print_ModelVariables(simulationModel::ModiaMath.AbstractSimulationModel) = print_ModelVariables(simulationModel.var)

getVariableName(model,vcat,vindex) = ModiaMath.DAE.getVariableName(model,vcat,vindex;
                                                                   xNames = model.var.x_names) 

getResultNames(model::Any) = model.var.result_names


function storeVariables!(res::ModiaMath.RawResult, model, t::Float64, x::Vector{Float64}, derx::Vector{Float64}, w::Vector{Float64})
    res.nt += 1
    i = res.nt
    ntMax = size(res.data,1)
    if i > ntMax
        # Result storage is full, double the current storage
        res.data = [res.data; zeros(ntMax,size(res.data,2))]
    end

    # Copy data in raw result storage
    for v in model.var.result_var
        if ModiaMath.isScalar(v)
            res.data[i, v.iresult] = v.value
        else
            vv = v.value
            for j in 1:length(vv)
                res.data[i, v.iresult+j-1 ] = vv[j]
            end
        end
    end
    return nothing
end

function getResult(model::ModiaMath.AbstractSimulationModel, res::ModiaMath.RawResult)
    nt         = res.nt
    data       = res.data

    seriesDict = Dict{AbstractString,Any}()
    varDict    = Dict{Symbol,ModiaMath.AbstractVariable}()
    for v in model.var.result_var
        name    = instanceName(v)
        len     = length(v.value)
        iresult = v.iresult
        if len == 1
            seriesDict[string(name)] = @view data[1:nt, iresult]
        else
            seriesDict[string(name)] = @view data[1:nt, iresult:iresult+len-1]
        end
        varDict[name] = v
   end
   result = ModiaMath.ResultWithVariables(seriesDict, varDict, model.modelName)
   return result, nt
end



computeVariables!(model::ModiaMath.AbstractComponentWithVariables,
                  simulationState::ModiaMath.SimulationState) = 
          error("Function ModiaMath.computeVariables!(model) not defined for model ", typeof(model))

function getModelResidues!(m::SimulationModel, t::Float64, _x::Vector{Float64}, _derx::Vector{Float64}, _r::Vector{Float64}, _w::Vector{Float64})  
    # Copy _x and _derx values to variables
    ModiaMath.copy_x_and_derx_to_variables!(t, _x, _derx, m.var)

    # Compute model variables
    m.simulationState.time = t
    computeVariables!(m.model, m.simulationState)

    # Copy variables to residues
    ModiaMath.copy_variables_to_residue!(m.var,_x,_derx,_r)  
    return nothing
end


# Return a table of actual variable and residue values from nonlinear solver in case of error
function ModiaMath.getVariableAndResidueValues(m::SimulationModel)
   var = m.var
   v_table = get_variableValueTable(var)
   r_table = get_residueValueTable(var, m.simulationState.derxev)
   return (v_table, r_table)
end
