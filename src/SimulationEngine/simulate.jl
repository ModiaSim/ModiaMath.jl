# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module
#   ModiaMath.SimulationEngine(ModiaMath/SimulationEngine/_module.jl)
#


mutable struct IntegratorData
    model::ModiaMath.AbstractSimulationModel    # Model specific data structure
    simulationState::ModiaMath.SimulationState
    statistics::ModiaMath.SimulationStatistics   # Simulation statistics
    z::Vector{Float64}          # Vector of event indicators (zero crossings). If one of z[i] passes
                                # zero, that is beforeEvent(z[i])*z[i] < 0, an event is triggered
                                # (allocated during instanciation according to nz).
    zDir::Vector{Cint}          # zDir[i] =  0: Root is reported for both crossing directions
                                #         =  1: Root is reported when crossing from negative to positive direction
                                #         = -1: Root is reported when crossing from positive to negative direction

    eweight::Sundials.N_Vector   # Weight vector
    last_t::Float64              # Time instant of the last call to idasol_f (residue function)
    last_norm_y::Float64         # Norm of y of the last call to idasol_f (residue function)
    hcur::MVector{1,Cdouble}     # Current step size argument for IDAGetCurrentStep
    order::MVector{1,Cint}       # Current order for IDAGetCurrentOrder
    r::Vector{Sundials.realtype} # Residues vector used in idasol_f function

    ida_mem::Ptr{Nothing}              # IDA pointer (to access all IDAgetXXX functions)

    y::Vector{Sundials.realtype}    # Julia vector of y wrapping y_N_Vector vector
    yp::Vector{Sundials.realtype}   # Julia vector of yp wrapping yp_N_Vector vector

    function IntegratorData(model::ModiaMath.AbstractSimulationModel, simulationState::ModiaMath.SimulationState)
        ny = simulationState.nx
        nz = simulationState.nz
        new(model, simulationState, simulationState.statistics,
          ones(nz), fill(0, 0), Sundials.N_VNew_Serial(ny), typemin(Float64), 0.0,
          MVector{1,Cdouble}(0.0), MVector{1,Cint}(0), zeros(ny))
    end
end

include("sundials_ida_UserData.jl")

function updateIDAstatistics!(mem, stat::ModiaMath.SimulationStatistics)
    info = Clong[0]
    Sundials.IDAGetNumSteps(mem, info)
    stat.nSteps        += info[1]

    stat.sparseSolver ? Sundials.IDASlsGetNumJacEvals(mem, info) : Sundials.IDADlsGetNumJacEvals(mem, info)
    stat.nJac          += info[1]

    Sundials.IDAGetNumErrTestFails(mem, info)
    stat.nErrTestFails += info[1]

    h0 = [0.0]
    Sundials.IDAGetActualInitStep(mem, h0)
    stat.h0 = min(stat.h0, h0[1])
end


"""
    result = closeTimePoints(t1,t2,epsilon)

The function returns true if the two time points t1 and t2 are close together (with respect to epsilon)
"""
closeTimePoints(t1::Float64, t2::Float64, epsilon::Float64) = abs(t1 - t2) / max(abs(t2), 10 * epsilon) <= 100 * epsilon


"""
    ModiaMath.simulate!(simulationModel; log=false, startTime=NaN, stopTime=NaN,
                                                    tolerance=NaN, interval=NaN)

Simulates a DAE `simulationModel` that is defined with package `Modia`, package `Modia3D` or
with the `ModiaMath.@component` macro. The DAE is mathematically described as
implicit, index 1 DAE with events (containing an ODE or a semi-explicit
index 1 DAE as special cases):

```math
\\begin{align}
 0       &= f_d(\\dot{x}, x, t, z_{pos}) \\\\
 0       &= f_c(x, t, z_{pos}) \\\\
 z       &= f_z(x,t) \\\\
 z_{pos} &= \\text{event}() \\; ? \\; z > 0 \\; : \\; \\text{last}(z_{pos}) \\\\
 J       &= \\left[ \\frac{\\partial f_d}{\\partial \\dot{x}};
                    \\frac{\\partial f_c}{\\partial x} \\right] \\; \\text{is regular}
\\end{align}
```

with initial conditions ``x_0^{-}``:

```math
\\lim_{\\epsilon \\rightarrow 0} x(t_0 - \\epsilon) = x_0^{-}
```

During continuous integration, equation system (1)-(4) is solved with the
[Sundials](https://computation.llnl.gov/projects/sundials) IDA solver
(accessed via the Julia [Sundials](https://github.com/JuliaDiffEq/Sundials.jl) interface package).
ModiaMath assumes that ``J`` (5) is **regular** for all time instants.
If this condition is violated, initialization and simulation will usually fail and an error message of the
form *"Solver does not converge"* might appear. Note, ModiaMath does not check this condition and can therefore
not provide better diagnostics in such cases.

If one of the elements of ``z`` crosses zero, an event is triggered and simulation is halted.
At an event, equation ``z_{pos} = z > 0`` (element wise) is added. The equation system (1)-(4)
is then solved with a fixed-point iteration scheme (= *event iteration*). Afterwards, integration is
restarted and ``z_{pos}`` keeps its value until the next event occurs.

Initial conditions ``x_0^{-}`` must be provided before simulation can start.
Best is if they fulfil the constraint equation ``0 = f_c(x_0^{-}, t_0, z > 0)``.
If this is not the case, initialization will simulate for an infinitesimal small time instant
so that ``x_0^{-}`` changes discontinuously to ``x_0^{+}`` with ``f_c (x_0^{+}, t_0, z > 0 )=0``.
Note, ``\\dot{x}`` is a Dirac impulse in this case.


Input arguments

- `simulationModel::ModiaMath.AbstractSimulationModel`: Model struct (generated with Modia, Modia3D or `ModiaMath.@component`).
- `log::Bool`: = true, if logging is enabled, otherwise it is disabled.
- `startTime::Float64`: Start time of the simulation in [s].
                        If startTime=NaN, the default startTime is used that is defined by the `simulationModel`.
- `stopTime::Float64`: Stop time of the simulation in [s].
                       If stopTime=NaN, the default stopTime is used that is defined by the `simulationModel`.
- `tolerance::Float64`: The relative tolerance for the integration.
  The absolute tolerance is computed as `0.1*tolerance*nominal(variable)` where
  `nominal(variable)` is the nominal value of the variable.
   If tolerance=NaN, the default tolerance is used that is defined by the `simulationModel`.
- `interval::Float64`: Output interval for results in [s]. If events occur, the event time instants are
  additionally added to the result.
  If interval=NaN, the default interval is used that is defined by the `simulationModel`.
"""
function simulate!(model::ModiaMath.AbstractSimulationModel;
                   tolerance=NaN,
                   startTime=NaN,
                   stopTime=NaN,
                   interval=NaN,
                   log::Bool=false,
                   KLUorderingChoice::Int=1)

    sim           = model.simulationState
    sim.model     = model
    sim.tolerance = isnan(tolerance) ? sim.defaultTolerance : convert(Float64, tolerance)
    sim.startTime = isnan(startTime) ? sim.defaultStartTime : convert(Float64, startTime)
    sim.stopTime  = isnan(stopTime)  ? sim.defaultStopTime  : convert(Float64, stopTime)
    sim.interval  = isnan(interval)  ? sim.defaultInterval  : convert(Float64, interval)
    sim.time      = sim.startTime
    @assert(sim.stopTime >= sim.startTime)
    @assert(sim.interval >= 0.0)
    @assert(sim.tolerance > 0.0)

    stopTime2::Float64  = sim.stopTime
    startTime2::Float64 = sim.startTime
    tolerance::Float64  = sim.tolerance
    interval::Float64   = sim.interval
    logger              = sim.logger
    ModiaMath.setLog!(logger, log)

    #if ModiaMath.isLogInfos(logger)
    println("... ModiaMath.simulate! (version 0.2.2-dev.2 from 2018-09-09 12:37) to simulate model: ", sim.name)
    #end

    # Start timing measure
    cpuStart::UInt64 = time_ns()
    cpuLast::UInt64  = cpuStart
    cpuStartIntegration::UInt64 = cpuStart

    # Init variables
    nt::Int           = round(Int64, 1 + (stopTime2 - startTime2) / interval)
    t0::Float64       = startTime2
    tEnd::Float64     = stopTime2
    tReached::Float64 = t0
    epsilon::Float64  = eps()

    # Initialize simulation model and store result after initialization
    statistics = sim.statistics
    ModiaMath.reInitializeStatistics!(statistics, t0, stopTime2, interval, tolerance)

    if ModiaMath.isLogInfos(logger)
        println("      Initialization at time = ", t0, " s")
    end

    init = DAE.initialize!(model, sim, t0, nt, tolerance)
    @assert(length(init.y0)  == sim.nx)
    @assert(length(init.yp0) == sim.nx)
    @assert(tolerance > 0.0)

    if init.terminate || tEnd <= t0
        statistics.h0   = 0.0
        statistics.hMin = 0.0
        @goto TerminateSimulation
    end

    if ModiaMath.isLogInfos(logger)
        println("      Simulation started")
    end

    maxTime::Float64          = init.maxTime
    nextEventTime::Float64    = init.nextEventTime
    integrateToEvent::Bool = init.integrateToEvent

    if nextEventTime < t0
        nextEventTime = typemax(Float64)
    end

    # Initialize auxiliary variables
    simModel = IntegratorData(model, sim)
    stateEvent::Bool = false
    timeEvent::Bool  = false
    flag::Cint = 0
    ny = sim.nx
    nz = sim.nz

    y  = copy(init.y0)
    yp = copy(init.yp0)
    simModel.y  = y
    simModel.yp = yp

    eventInfo = DAE.EventInfo()
    tret    = [0.0]
    jc::Int = 0  # communication point counter
    if ModiaMath.isLogInfos(logger) && nz > 0
        rootInfo = fill(Cint(0), nz)
    end

    # Create IDA
    mem = Sundials.IDACreate()
    if mem == C_NULL
        error("ModiaMath.simulate: Failed to allocate IDA solver object")
    end
    simModel.ida_mem = mem

    # Allocate N_Vector storage for y and yp
    y_N_Vector  = Sundials.N_VMake_Serial(ny, pointer(y))
    yp_N_Vector = Sundials.N_VMake_Serial(ny, pointer(yp))

    # Run IDA
    try
       # Set error handler function
        Sundials.IDASetErrHandlerFn(mem, @cfunction(idasol_ErrHandlerFn, Nothing, (Cint,Cstring,Cstring,Cstring,Ref{typeof(IntegratorData)})), pointer_from_objref(simModel))

        # Initialize IDA
        Sundials.IDAInit(mem, @cfunction(idasol_f, Cint, (Sundials.realtype,Sundials.N_Vector,Sundials.N_Vector,Sundials.N_Vector,Ref{typeof(IntegratorData)})), t0, y_N_Vector, yp_N_Vector)

        #IDASStolerances(mem, tolerance, 0.1*tolerance)
        tolAbs = 0.1 * tolerance * init.y_nominal
        for i in 1:ny
            if !init.y_errorControl[i]
                tolAbs[i] = 1e5 * init.y_nominal[i]   # switch tolerance control off
            end
        end

        Sundials.IDASVtolerances(mem, tolerance, tolAbs)
        Sundials.IDASetUserData(mem, simModel)

        # Sundials.IDASetId(mem, y_states_nvector)
        #if sim.sparse
        #   nnz_jac = nnz(sim.jac)
        #   Sundials.IDAKLU(mem, ny, nnz_jac)
        #   Sundials.IDAKLUSetOrdering(mem, KLUorderingChoice)
        #   Sundials.IDASlsSetSparseJacFn(mem);
        #else
        A  = Sundials.SUNDenseMatrix(ny, ny)
        LS = Sundials.SUNDenseLinearSolver(y, A)
        Sundials.IDADlsSetLinearSolver(mem, LS, A)
        #end

        # Initialize zero crossing function, if required
        hasZeroCross = nz > 0
        if hasZeroCross
            Sundials.IDARootInit(mem, nz, @cfunction(idasol_g, Cint, (Sundials.realtype,Sundials.N_Vector,Sundials.N_Vector,Ptr{Sundials.realtype},Ref{typeof(IntegratorData)})))
            Sundials.IDASetRootDirection(mem, sim.zDir)
        end

        # Initialize event variables
        tStop::Float64 = tEnd
        tNext = t0
        tc    = t0
        nTimeEvents = 0
        restart     = ModiaMath.NoRestart

        # Check for event at time = t0. If yes, trigger event and store result after event
        if closeTimePoints(t0, nextEventTime, epsilon)
            tReached = t0
            if ModiaMath.isLogEvents(logger)
                println("\n      Time event at time = $tReached s")
            end

            DAE.reset!(eventInfo)
            DAE.processEvent!(model, sim, tReached, y, yp, eventInfo)
            restart          = eventInfo.restart
            maxTime          = eventInfo.maxTime
            nextEventTime    = eventInfo.nextEventTime
            integrateToEvent = eventInfo.integrateToEvent

            statistics.nTimeEvents += 1
            if ModiaMath.isLogEvents(logger)
                println("        restart = ", restart)
            end

            # Reinitialize IDA if required by model
            if restart == ModiaMath.Restart
                statistics.nRestartEvents += 1
                Sundials.__IDAReInit(mem, tReached, y_N_Vector, yp_N_Vector);
            end
        end

        # Loop over communication points
        cpuStartIntegration = time_ns()

        while tReached < tEnd && restart != ModiaMath.Terminate  #------------------ while --------------
            if tReached >= tc
                # Last communication point reached -> increment communication point
                jc += 1
                tc = startTime2 + jc * interval
            end

            if closeTimePoints(nextEventTime, tc, epsilon)
                nextEventTime = tc
            end

            if nextEventTime <= tc # next tStop is a time event
                event = true
                tNext = nextEventTime
                tStop = integrateToEvent ? nextEventTime : min(tEnd, maxTime)
            else        # next tStop is a communication point (no event)
                event = false
                tNext = tc
                tStop = min(tEnd, maxTime)
            end

            # Integrate in direction of tStop
            #println("... tReached = ", tReached, ", tNext = ", tNext, ", tStop = ", tStop,
            #        ", y = ", y[1], ", yp = ", yp[1])
            Sundials.IDASetStopTime(mem, tStop)
            #flag = Sundials.IDASolve(mem, tNext, tret, simModel.y, simModel.yp, Sundials.IDA_NORMAL)
            flag = Sundials.__IDASolve(mem, tNext, tret, y_N_Vector, yp_N_Vector, Sundials.IDA_NORMAL)
            tReached = tret[1]
            #println("    tReached = ", tReached)

            stateEvent = flag == Sundials.IDA_ROOT_RETURN
            timeEvent  = event && closeTimePoints(tNext, tReached, epsilon)
            isEvent    = timeEvent || stateEvent

            # Store result point at tReached
            DAE.computeAndStoreResult!(model, sim, tReached, y, yp)

            # If event at tReached, handle the event
            if isEvent
                if ModiaMath.isLogEvents(logger)
                    if timeEvent && stateEvent
                        print("\n      Time and state (zero-crossing) event at time = $tReached s")
                    elseif timeEvent
                        println("\n      Time event at time = $tReached s")
                    elseif stateEvent
                        print("\n      State event (zero-crossing) at time = $tReached s")
                    end

                    if stateEvent
                       # Print information about the root
                        Sundials.IDAGetRootInfo(mem, rootInfo)
                        print(" (")
                        firstRoot = true
                        for iroot in eachindex(rootInfo)
                            if rootInfo[iroot] != 0
                                if firstRoot
                                    firstRoot = false
                                else
                                    print(", ")
                                end
                                print("z[", iroot, "]", rootInfo[iroot] > 0 ? " > 0" : " < 0")
                            end
                        end
                        print(")\n")
                    end
                end

                # Trigger event
                DAE.reset!(eventInfo)
                DAE.processEvent!(model, sim, tReached, y, yp, eventInfo)
                restart          = eventInfo.restart
                maxTime          = eventInfo.maxTime
                nextEventTime    = eventInfo.nextEventTime
                integrateToEvent = eventInfo.integrateToEvent

                if timeEvent
                    statistics.nTimeEvents += 1
                end

                if stateEvent
                    statistics.nStateEvents += 1
                end

                if restart == ModiaMath.Restart || restart == ModiaMath.FullRestart
                    statistics.nRestartEvents += 1
                end

                if ModiaMath.isLogEvents(logger)
                    println("        restart = ", restart)
                end

                # Handle restart flag
                if restart == ModiaMath.Restart
                    # Restart (dimensions do not change)
                    updateIDAstatistics!(mem, statistics)
                    Sundials.__IDAReInit(mem, tReached, y_N_Vector, yp_N_Vector);
                elseif restart == ModiaMath.FullRestart
                    error("FullRestart not yet supported")
                end
            end

            # print every 5 s the time in order that the user sees the progress
            if ModiaMath.isLogProgress(logger)
                cpuNew = time_ns()
                if (cpuNew - cpuLast) * 1e-9 > 5.0
                    cpuLast = cpuNew
                    @printf("      progress: integrated up to time = %.2g s\n", tReached)
                end
            end
        end #---------------------- end while ------------------------------------

        # Finalize statistics for simulation run
        updateIDAstatistics!(mem, statistics)
    finally
        # Free allocated memory
        Sundials.__N_VDestroy_Serial(y_N_Vector)
        Sundials.__N_VDestroy_Serial(yp_N_Vector)
        Sundials.IDAFree([mem])
    end

    @label TerminateSimulation
    if ModiaMath.isLogInfos(logger)
        println("\n      Simulation is terminated at time = ", tEnd, " s")
    end

    result = DAE.terminate!(model, sim, tReached, y, yp)
    statistics.cpuTimeInitialization = (cpuStartIntegration - cpuStart) * 1e-9
    statistics.cpuTimeIntegration    = (time_ns() - cpuStartIntegration) * 1e-9

    if ModiaMath.isLogStatistics(logger)
        println("\n      Statistics (get help with ?ModiaMath.SimulationStatistics):")
        display(statistics)
    end
    return result
end
