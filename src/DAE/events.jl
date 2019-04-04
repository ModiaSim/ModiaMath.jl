# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.DAE (ModiaMath/DAE/_module.jl)
#

# Provide functions to handle time and state events in user models.
@eval using Printf

mutable struct EventHandler
   # Input values for the event functions
    time::Float64              # Current simulation time
    initial::Bool              # = true, if model is called at initialization
                              #         (if initial, event=true)
    terminal::Bool             # = true, if model is called for termination (close files, streams, visualization, ...)
    event::Bool                # = true, if model is called at an event  
    afterSimulationStart::Bool # = true, if model is called after simulation start
    crossing::Bool             # = true, if model is called to compute crossing function
   
    # Computed by the event functions
    # For time events:
    maxTime::Float64        # Integrate at most up to maxTime
                            # (if integrateToEvent == true, maxTime = nextEventTime,
                            #  otherwise maxTime > nextEventTime)
    nextEventTime::Float64  # Next time event instant (= realmax(Time) if no event)
    integrateToEvent::Bool  # = true, if integrator shall integrate to nextEventTime
                            # = false, if integrator can integrate beyond nextEventTime 
    restart::EventRestart   # Type of integrator restart at current event
    newEventIteration::Bool   # = true, if another event iteration; = false, if no event iteration anymore
    firstEventIteration::Bool # = true, if first iteration at an event.
   
   # For state events:
    nz::Int                 # Number of event indicators
    z::Vector{Float64}      # Vector of event indicators (zero crossings). If one of z[i] passes
                            # zero, that is beforeEvent(z[i])*z[i] < 0, an event is triggered
                            # (allocated during instanciation according to nz).
    zPositive::Vector{Bool} # = true if z > 0 at the last event instant, otherwise false
    zDir::Vector{Cint}      # zDir[i] =  0: Root is reported for both crossing directions
                            #         =  1: Root is reported when crossing from negative to positive direction
                            #         = -1: Root is reported when crossing from positive to negative direction
   
    # Logger
    logger::ModiaMath.Logger
     
    function EventHandler(;nz::Int=0) 
        @assert(nz >= 0)
        new(0.0, false, false, false, false, false, floatmax(Float64), floatmax(Float64),
            true, NoRestart, false, false, nz, ones(nz), fill(false, nz), fill(0, nz))
    end
end

positiveCrossingAsString(positive::Bool) = positive ? " (became > 0)" : " (became <= 0)"
negativeCrossingAsString(negative::Bool) = negative ? " (became < 0)" : " (became >= 0)"

   
function initEventIteration!(h::EventHandler, t::Float64)
    h.time          = t
    h.restart       = NoRestart
    h.maxTime       = floatmax(Float64)
    h.nextEventTime = floatmax(Float64)
    h.newEventIteration   = false
    h.firstEventIteration = true
end

function terminateEventIteration!(h::EventHandler)
    result                = h.initial ? false : !h.newEventIteration
    h.firstEventIteration = false
    h.newEventIteration   = false
    return result
end


isInitial(h::EventHandler)              = h.initial
isTerminal(h::EventHandler)             = h.terminal
isEvent(h::EventHandler)                = h.event
isAfterSimulationStart(h::EventHandler) = h.afterSimulationStart
isEventIndicator(h::EventHandler) = h.crossing
const zEps = 1.e-14


function setNextEvent!(h::EventHandler, nextEventTime::Float64; 
                       integrateToEvent::Bool=true, 
                       restart::EventRestart=Restart)
    if (h.event && nextEventTime > h.time) || (h.initial && nextEventTime >= h.time)
        if integrateToEvent
            h.maxTime = min(h.maxTime, nextEventTime)
        end
        if nextEventTime < h.nextEventTime
            h.nextEventTime    = nextEventTime
            h.integrateToEvent = integrateToEvent
            if ModiaMath.isLogEvents(h.logger)
                @printf("        nextEventTime = %.6g s, integrateToEvent = %s\n", 
                   nextEventTime, integrateToEvent ? "true" : "false")
            end
        end
        h.restart = max(h.restart, restart)
    end
    return nothing
end


function positive!(h::EventHandler, nr::Int, crossing::Float64, crossingAsString::String; 
                   restart::EventRestart=Restart)::Bool
    if h.initial
        h.zPositive[nr] = crossing > 0.0
        if ModiaMath.isLogEvents(h.logger)
            println("        ", crossingAsString, " = ", crossing, positiveCrossingAsString(h.zPositive[nr]))
        end

    elseif h.event
        # println("... nr = ", nr, ", crossing = ", crossing)
        new_zPositive = crossing > 0.0
        change = (h.zPositive[nr] && !new_zPositive) || (!h.zPositive[nr] && new_zPositive)
        h.zPositive[nr] = new_zPositive
        
        if change
            h.restart = max(h.restart, restart)
            if ModiaMath.isLogEvents(h.logger)
                println("        ", crossingAsString, " = ", crossing, positiveCrossingAsString(h.zPositive[nr]))
            end
            h.newEventIteration = true
        end
    end
    h.z[nr] = crossing + (h.zPositive[nr] ? zEps : -zEps)
   
    return h.zPositive[nr]
end


function negative!(h::EventHandler, nr::Int, crossing::Float64, crossingAsString::String; 
                   restart::EventRestart=Restart)::Bool
    if h.initial
        h.zPositive[nr] = crossing >= 0.0
        if ModiaMath.isLogEvents(h.logger)
            println("        ", crossingAsString, " = ", crossing, negativeCrossingAsString(!h.zPositive[nr]))
        end

    elseif h.event
        new_zPositive = crossing >= 0.0
        change = (h.zPositive[nr] && !new_zPositive) || (!h.zPositive[nr] && new_zPositive)
        h.zPositive[nr] = new_zPositive
        
        if change
            h.restart = max(h.restart, restart)
            if ModiaMath.isLogEvents(h.logger)
                println("        ", crossingAsString, " = ", crossing, negativeCrossingAsString(!h.zPositive[nr]))
            end
            h.newEventIteration = true
        end
    end
    h.z[nr] = crossing + (h.zPositive[nr] ? zEps : -zEps)
  
    return !h.zPositive[nr]
end


function change!(h::EventHandler, nr::Int, crossing::Float64, crossingAsString::String; 
                 restart::EventRestart=Restart)::Bool
    h.z[nr] = crossing
    if h.initial 
        h.zPositive[nr] = crossing > 0.0
        if ModiaMath.isLogEvents(h.logger)
            println("        ", crossingAsString, " = ", crossing, positiveCrossingAsString(h.zPositive[nr]))
        end      
        return false

    elseif h.event
        new_zPositive = crossing > 0.0
        change = (h.zPositive[nr] && !new_zPositive) || (!h.zPositive[nr] && new_zPositive)
        h.zPositive[nr] = new_zPositive
        
        if change
            h.restart = max(h.restart, restart)
            if ModiaMath.isLogEvents(h.logger)
                println("        ", crossingAsString, " = ", crossing, positiveCrossingAsString(h.zPositive[nr]))
            end
            h.newEventIteration = true        
            return true
        end
    end  

    return false
end


function edge!(h::EventHandler, nr::Int, crossing::Float64, crossingAsString::String; 
               restart::EventRestart=Restart)::Bool
    h.z[nr] = crossing
    
    if h.initial 
        h.zPositive[nr] = crossing > 0.0
        h.zDir[nr] = 1
    
        if ModiaMath.isLogEvents(h.logger)
            println("        ", crossingAsString, " = ", crossing, positiveCrossingAsString(h.zPositive[nr]))
        end        
        return false
    
    elseif h.event
        new_zPositive = crossing > 0.0
        edge = !h.zPositive[nr] && new_zPositive
        h.zPositive[nr] = new_zPositive
    
        if edge
            h.restart = max(h.restart, restart)
            if ModiaMath.isLogEvents(h.logger)
                println("        ", crossingAsString, " = ", crossing, " (became > 0)")
            end
            h.newEventIteration = true        
            return true
        end

    elseif h.crossing && h.zPositive[nr] && crossing < 0.0
        # Eventually need to improve later
        # (need to check, whether h.get_z is good enough to switch h.zPositive)
        h.zPositive[nr] = false
    end
    return false
end
