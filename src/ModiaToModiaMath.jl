# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

module ModiaToModiaMath

import ModiaMath
using ModiaMath.DAE

mutable struct ModiaSimulationModel <: ModiaMath.AbstractSimulationModel
    modelName::String
    simulationState::ModiaMath.SimulationState
    getModiaResidues!::Function
    preInitial::Bool
    nz_preInitial::Int
    nw::Int 
    store    # extra storage
   
    function ModiaSimulationModel(name::Symbol, f!::Function, x0::Vector{Float64}, der_x0::Vector{Float64}, jac=nothing;
                                structureOfDAE = ModiaMath.DAE_NoSpecialStructure,
                                is_constraint=fill(false,length(x0)),
								has_constraintDerivatives=false,
                                xNames::Vector{String}=ModiaMath.DAE.nameVector("x", length(x0)),
                                derxNames::Vector{String}=ModiaMath.DAE.fcNameVector("der", xNames),
                                wNames::Vector{String}=fill("", 0),
                                xw_states=fill(true, length(x0)),
                                nc=0, nz=0, maxSparsity=0.1,
                                startTime=0.0,
                                stopTime=1.0,
                                tolerance=1e-4,
                                interval=(stopTime - startTime) / 500.0)                

        simulationState = ModiaMath.SimulationState(string(name), getModelResidues!, x0;
                                structureOfDAE = structureOfDAE,
                                is_constraint = is_constraint,
								has_constraintDerivatives = has_constraintDerivatives,
                                nc=nc == 0 ? 1 : nc,
                                nz=nz,
                                jac=jac,
                                scaleConstraintsAtEvents=false,
                                maxSparsity=maxSparsity,
                                defaultStartTime=startTime,
                                defaultStopTime=stopTime,
                                defaultTolerance=tolerance,
                                defaultInterval=interval)
                                
        new(string(name), simulationState, f!, false, 0, length(wNames), nothing)     
    end      

    function ModiaSimulationModel(name::String, 
                            getModelResidues!::Function,
                            x_start::Vector{Float64},
                            getVariableName::Function=ModiaMath.DAE.defaultVariableName; 
                            structureOfDAE = ModiaMath.DAE_NoSpecialStructure,     
                            is_constraint=fill(false,length(x_start)),
							has_constraintDerivatives=false,							
                            nc::Int=1,
                            nz::Int=0,                   
                            nw::Int=0,
                            x_fixed::Vector{Bool}=fill(false, length(x_start)),      
                            x_nominal::Vector{Float64}=fill(NaN, length(x_start)),                                                 
                            hev::Float64=1e-8,
                            scaleConstraintsAtEvents::Bool=false,               
                            jac=nothing, maxSparsity::Float64=0.1,
                            startTime=0.0,
                            stopTime=1.0,
                            tolerance=1e-4,
                            interval=(stopTime - startTime) / 500.0)
      
        simulationState = ModiaMath.SimulationState(name, getModelResidues!, x_start, getVariableName;
                            structureOfDAE = structureOfDAE,
                            is_constraint = is_constraint, 
							has_constraintDerivatives = has_constraintDerivatives,
                            nc=nc,
                            nz=nz,
                            nw=nw,
                            x_fixed=x_fixed,
                            x_nominal=x_nominal,
                            hev=hev,
                            scaleConstraintsAtEvents=scaleConstraintsAtEvents,
                            jac=jac,
                            maxSparsity=maxSparsity,
                            defaultStartTime=startTime,
                            defaultStopTime=stopTime,
                            defaultTolerance=tolerance,
                            defaultInterval=interval)
                                
        new(name, simulationState, getModelResidues!, false, 0, nw, nothing) 
    end
   

    function ModiaSimulationModel()
        function emptyFunction
        end

        simulationState = ModiaMath.DAE.SimulationState("???", emptyFunction, zeros(1))
        new("???", simulationState, emptyFunction, true, 0, 0, nothing)
    end
end 


function getModelResidues!(m::ModiaSimulationModel, t::Float64, y::Vector{Float64}, yp::Vector{Float64}, residues::Vector{Float64}, w::Vector{Float64})
    m.time = t
    if m.nw == 0
        Base.invokelatest(m.getModiaResidues!, m, y, yp, residues)
    else
        Base.invokelatest(m.getModiaResidues!, m, y, yp, residues, w)
    end 	  
end


function getRawResult(m::ModiaSimulationModel)
    sim   = m.simulationState
    res   = sim.rawResult
    nt    = res.nt
    nx    = sim.nx
    t     = @view res.data[1:nt,1]
    X     = @view res.data[1:nt,2:nx + 1]
    DER_X = @view res.data[1:nt,nx + 2:2 * nx + 1]
    W     = @view res.data[1:nt,2 * nx + 2:end]
    return (t, X, DER_X, W)
end
 

#export ModiaSimulationModel
#export EventRestart, NoRestart, Restart, FullRestart, Terminate
#export isPreInitial, isInitial, isEvent, isLog
#export setNextEvent!, positive!, positiveChange!, positiveEdge!

isPreInitial(m::ModiaSimulationModel) = m.preInitial
isInitial(m::ModiaSimulationModel) = ModiaMath.isInitial(m.simulationState)
isTerminal(m::ModiaSimulationModel) = ModiaMath.isTerminal(m.simulationState)
isStoreResult(m::ModiaSimulationModel) = ModiaMath.isStoreResult(m.simulationState)
isLogInfos(m::ModiaSimulationModel) = ModiaMath.isLogInfos(m.simulationState)
isLogWarnings(m::ModiaSimulationModel) = ModiaMath.isLogWarning(m.simulationState)
isLogEvents(m::ModiaSimulationModel) = ModiaMath.isLogEvents(m.simulationState)
isEvent(m::ModiaSimulationModel) = m.simulationState.eventHandler.event
isAfterSimulationStart(m::ModiaSimulationModel) = m.simulationState.eventHandler.afterSimulationStart


setNextEvent!(m::ModiaSimulationModel, nextEventTime::Float64;
              integrateToEvent::Bool=true, restart::ModiaMath.EventRestart=ModiaMath.Restart) = ModiaMath.setNextEvent!(m.simulationState.eventHandler, nextEventTime; 
                                                                                                     integrateToEvent=integrateToEvent, restart=restart)
    
function positive!(crossing::Float64, m::ModiaSimulationModel, nr::Int, crossingAsString="????"; restart::ModiaMath.EventRestart=ModiaMath.Restart)
    if m.preInitial
        m.nz_preInitial = max(m.nz_preInitial, nr)
        return crossing > 0.0
    end
    ModiaMath.DAE.positive!(m.simulationState.eventHandler, nr, crossing, crossingAsString; restart=restart)
end

function positiveChange!(crossing::Float64, m::ModiaSimulationModel, nr::Int, crossingAsString="????"; restart::ModiaMath.EventRestart=ModiaMath.Restart)
    if m.preInitial
        m.nz_preInitial = max(m.nz_preInitial, nr)
        return false
    end
    ModiaMath.DAE.change!(m.simulationState.eventHandler, nr, crossing, crossingAsString="????"; restart=ModiaMath.Restart)
end

function positiveEdge!(crossing::Float64, m::ModiaSimulationModel, nr::Int, crossingAsString="????"; restart::ModiaMath.EventRestart=ModiaMath.Restart)
    if m.preInitial
        m.nz_preInitial = max(m.nz_preInitial, nr)
        return false
    end
    ModiaMath.DAE.edge!(m.simulationState.eventHandler, nr, crossing, crossingAsString="????"; restart=ModiaMath.Restart)
end


simulate(m::ModiaSimulationModel, t; log::Bool=false, tolRel::Float64=1e-4, KLUorderingChoice::Int=1) =
           simulate(m, collect(t); log=log, tolRel=tolRel, KLUorderingChoice=KLUorderingChoice)

function simulate(m::ModiaSimulationModel, t::Vector{Float64};
                  log::Bool=false, 
                  tolRel::Float64=1e-4, 
                  KLUorderingChoice::Int=1)     
    startTime = t[1]
    stopTime = t[end]   
    interval = (stopTime - startTime) / (length(t) - 1)  
    
    ModiaMath.simulate!(m;stopTime=stopTime, 
                        tolerance=tolRel, 
                        startTime=startTime, 
                        interval=interval, 
                        log=log,
                        KLUorderingChoice=KLUorderingChoice)
end                    
                   
end