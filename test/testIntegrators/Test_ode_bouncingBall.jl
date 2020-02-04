module Test_ode_bouncingBall

#=
Try EventHandler of ModiaMath for DifferentialEquations.jl
to figure out what needs to be adapted.

Evaluate Measurements.jl for events (does not work yet in a good way; result is useless)

=#


using DifferentialEquations, Measurements, StaticArrays, PyPlot

MeanValue(v::Number) = v
MeanValue(v::Measurement{Float64}) = Measurements.value(v)

mutable struct EventHandler{M<:Real, F<:Real}   # M: Real type of variables where Measurement can be used.
                                                # F: Floating-point type for variables that are no measurements

   # Input values for the event functions
    time::M                    # Current simulation time
    initial::Bool              # = true, if model is called at initialization
                               #         (if initial, event=true)
    terminal::Bool             # = true, if model is called for termination (close files, streams, visualization, ...)
    event::Bool                # = true, if model is called at an event
    afterSimulationStart::Bool # = true, if model is called after simulation start
    crossing::Bool             # = true, if model is called to compute crossing function

    # Computed by the event functions
    # For time events:
    maxTime::F              # Integrate at most up to maxTime
                            # (if integrateToEvent == true, maxTime = nextEventTime,
                            #  otherwise maxTime > nextEventTime)
    nextEventTime::M        # Next time event instant (= typemax(Time) if no event)
    integrateToEvent::Bool  # = true, if integrator shall integrate to nextEventTime
                            # = false, if integrator can integrate beyond nextEventTime
    newEventIteration::Bool   # = true, if another event iteration; = false, if no event iteration anymore
    firstEventIteration::Bool # = true, if first iteration at an event.

   # For state events:
    nz::Int                 # Number of event indicators
    z::Vector{M}            # Vector of event indicators (zero crossings). If one of z[i] passes
                            # zero, that is beforeEvent(z[i])*z[i] < 0, an event is triggered
                            # (allocated during instanciation according to nz).
    zPositive::Vector{Bool} # = true if z > 0 at the last event instant, otherwise false
    zDir::Vector{Cint}      # zDir[i] =  0: Root is reported for both crossing directions
                            #         =  1: Root is reported when crossing from negative to positive direction
                            #         = -1: Root is reported when crossing from positive to negative direction

    function EventHandler{M,F}(;nz::Int=0) where {M<:Real, F<:Real}
        @assert(nz >= 0)
        new(0, false, false, false, false, false, M(floatmax(F)), M(floatmax(F)),
            true, false, false, nz, ones(M,nz), fill(false, nz), fill(0, nz))
    end
end


const zEps = 1e-12

function edge!(h::EventHandler, nr::Int, crossing)::Bool
    returnValue = false

    if h.initial
        h.zPositive[nr] = crossing > 0.0
        h.zDir[nr] = 1
        returnValue = false

    elseif h.event
        new_zPositive = crossing > 0.0
        edge = !h.zPositive[nr] && new_zPositive
        h.zPositive[nr] = new_zPositive

        if edge
            h.newEventIteration = true
            returnValue = true
        end

    elseif h.crossing && h.zPositive[nr] && crossing < 0.0
        # Eventually need to improve later
        # (need to check, whether h.get_z is good enough to switch h.zPositive)
        h.zPositive[nr] = false
    end
    h.z[nr] = crossing + (h.zPositive[nr] ? zEps : -zEps)
    return returnValue
end


mutable struct Model{M<:Real, F<:Real}   # M: Real type of variables where Measurement can be used.
                                         # F: Floating-point type for variables that are no measurements
    # parameters
    m::M   # mass
    e::M   # coefficient of restitution
    g::F   # gravity (not used as Measurment)

	# variables
	h::M  # height
	v::M  # velocity
	a::M  # acceleration
	flying::Bool   # = true, if flying

    # simulation state
	x_start::Vector{M}
	eventHandler::EventHandler

    function Model{M,F}(m::M, e::M, g::F, h0::M) where {M<:Real, F<:Real}
	    @assert(m >= 0)
	    @assert(e >= 0)
		@assert(h0 > 0)
		x_start = M[h0, 0]
        new(m, e, g, 0, 0, 0, true, x_start, EventHandler{M,F}(nz=1))
    end
end

Model(;m::M=1.0, e::M=0.7, g::F=9.81, h0::M=1.0) where {M<:Real, F<:Real} = Model{M,F}(m,e,g,h0)


"Compute m::Model at actual time instant"
function equations!(m::Model{M,F}, t, x)::Nothing where {M,F}
    m.h = x[1]
    m.v = x[2]
	if edge!(m.eventHandler, 1, -m.h)
		m.v = -m.e*m.v
		if m.v < 0.01
			m.flying = false
			m.v = 0.0
		end
		x[2] = m.v   # re-initialize state vector x
	end
	m.a = m.flying ? -m.g : 0.0
    return nothing
end

"derivatives!: Called by integrator"
function derivatives!(der_x, x, m::Model, t)::Nothing
    equations!(m, t, x)
    der_x[1] = m.v
    der_x[2] = m.a
    return nothing
end

"outputs!: Called at every communication point"
function outputs!(x::Vector{M}, t, integrator)::Model where {M<:Real}
    m::Model = integrator.p
    equations!(m, t, x)
    return deepcopy(m)
end


function outputs2!(x, t, integrator)::Nothing
    # println("outputs!2 called at time = ", t)
    return nothing
end

"condition_neg!: Called to compute the zero crossings from positive to negative"
function conditions_neg!(z, x, t, integrator)::Nothing
    m::Model = integrator.p
    h = m.eventHandler
    h.crossing = true
    equations!(m, t, x)
    h.crossing = false
    z[1] = -h.z[1]
    #println("conditions! called at time = ", t, ", z = ", z[1])
    return nothing
end

"affect!: Called at an event instant (either as affect! or here as affect_neg!)"
function affect!(integrator, event_index)::Nothing
    println("affect! called at time = ", integrator.t, ", h = ", integrator.u[1])
    m::Model = integrator.p
    h = m.eventHandler
    h.event = true
    equations!(m, integrator.t, integrator.u)
    h.event = false
	set_proposed_dt!(integrator, 1e-5)
    return nothing
end

function x₀(m::Model{M,F}, t0::M)::Vector{M} where {M<:Real, F<:Real}
    # Initialize model
    return deepcopy(m.x_start)
end


t_start = 0.0
t_inc   = 0.0001
t_end   = 5.0
tspan   = (t_start, t_end)
tspan2  = t_start:t_inc:t_end


# Problem 1 with Float64
model1 = Model()
println("typeof(model1) = ", typeof(model1))

# x1₀ = deepcopy(model1.x_start)

saved_values1=SavedValues(typeof(model1.x_start[1]), Model)
cb1 = SavingCallback(outputs!, saved_values1, saveat=tspan2)
cb2 = VectorContinuousCallback(conditions_neg!, nothing, 1; affect_neg! = affect!, interp_points=10)
cb3 = FunctionCallingCallback(outputs2!, funcat=tspan2)

prob1 = ODEProblem(derivatives!, x₀, tspan, model1)
sol1  = solve(prob1, Tsit5(), reltol = 1e-6, saveat=t_inc, callback=CallbackSet(cb1,cb2,cb3))

println("Simulation successful")

t = saved_values1.t
m = saved_values1.saveval
flying = getfield.(m, :flying)

figure(1)
clf()
plot(sol1.t, getindex.(sol1.u, 1), label="h")
plot(t, flying, label="flying")
grid(true)
legend()
title("Bouncing ball with Float64")

println("length(sol1.t) = ", length(sol1.t))
println("length(t) = ", length(t))


# Problem 2 with Measurement
Float64(v::Measurement{Float64}) = Measurements.value(v)

model2 = Model(m=1.0 ± 0.1, e=0.7 ± 0.1, h0=1.0± 0.1, g=9.81)
#println("model2 = $model2")
#println("typeof(model2) = ", typeof(model2))
x2₀ = deepcopy(model2.x_start)

t_start = 0.0 ± 0.0
t_inc   = 0.0001 ± 0.0
t_end   = 5.0 ± 0.0
tspan   = (t_start, t_end)
tspan2  = t_start:t_inc:t_end

saved_values2=SavedValues(typeof(model2.x_start[1]), Model)
cb21 = SavingCallback(outputs!, saved_values2, saveat=tspan2)
cb22 = VectorContinuousCallback(conditions_neg!, nothing, 1; affect_neg! = affect!, interp_points=10)
prob2 = ODEProblem(derivatives!, x2₀, tspan, model2)
sol2  = solve(prob2, Tsit5(), reltol = 1e-6, saveat=t_inc, callback=CallbackSet(cb21,cb22))

h = getindex.(sol2.u, 1)
hm = Measurements.value.(h)
hu = Measurements.uncertainty.(h)
hmax = hm + hu
hmin = hm - hu

figure(2)
clf()
plot(sol2.t, hmax, label="\$h_{max}\$")
plot(sol2.t, hmin, label="\$h_{min}\$")
plot(sol2.t, h, label="\$h_{mean}\$")
grid(true)
legend()
title("Bouncing ball with Measurement{Float64}")
#@show sol2.t
end