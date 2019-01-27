# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.DAE (ModiaMath/DAE/_module.jl)
#


"""
    initInfo = InitInfo(...)

Data structure that is returned when calling initialize!(...)
"""
struct InitInfo
    y0::Vector{Float64}          # Initial y
    yp0::Vector{Float64}         # Initial yp0
    y_nominal::Vector{Float64}   # Nominal values of y (y_nominal[i]>0.0 required)
    y_errorControl::Vector{Bool} # = true: error control; = false, no error control
    maxTime::Float64             # Integrate at most until maxTime
    nextEventTime::Float64       # Next time event
    integrateToEvent::Bool       # = true, if exact integration to nextEventTime
    terminate::Bool              # = true, terminate simulation

    function InitInfo(y0::Vector{Float64},
                      yp0::Vector{Float64};
                      y_nominal::Vector{Float64}=ones(length(y0)),
                      y_errorControl::Vector{Bool}=fill(true, length(y0)),                   
                      maxTime::Float64=floatmax(Float64),
                      nextEventTime::Float64=floatmax(Float64),
                      integrateToEvent::Bool=true,
                      terminate::Bool=false)
        @assert(length(y0) == length(yp0))
        @assert(length(y0) == length(y_nominal))
        @assert(y_nominal[ argmin(y_nominal) ] > 0.0)    
 
        new(y0, yp0, y_nominal, y_errorControl, maxTime, nextEventTime, integrateToEvent, terminate)
    end
end


"""
    eventInfo = EventInfo(...)

Data structure in which the return values of processEvent!(...) are stored
"""
mutable struct EventInfo
    restart::EventRestart
    maxTime::Float64         # Integrate at most until maxTime
    nextEventTime::Float64   # Next time event
    integrateToEvent::Bool   # = true, if exact integration to nextEventTime
    nextModel                # If restart==FullRestart, next model

    EventInfo() = new(Restart, floatmax(Float64), floatmax(Float64), true, nothing)
end


function reset!(einfo::EventInfo)
    einfo.restart = Restart
    einfo.maxTime = floatmax(Float64)
    einfo.nextEventTime = floatmax(Float64)
    einfo.integrateToEvent = true
    einfo.nextModel = nothing
end

