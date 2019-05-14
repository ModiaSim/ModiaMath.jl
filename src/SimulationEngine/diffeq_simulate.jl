# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module
#   ModiaMath.SimulationEngine(ModiaMath/SimulationEngine/_module.jl)
#
import DASSL, Sundials
using DASSL, Sundials

mutable struct SolData
       abstol::Float64
       reltol::Vector{Float64}
       hcur::Float64
       hmin::Float64
       hmax::Float64
       order::Int32
       weights::Vector{Float64}
       steps::Int32
       jac_evals::Int32
       fails::Int32
       end

mutable struct IntegratorData
    model::ModiaMath.AbstractSimulationModel    # Model specific data structure
    simulationState::ModiaMath.SimulationState
    statistics::ModiaMath.SimulationStatistics   # Simulation statistics
    z::Vector{Float64}          # Vector of event indicators (zero crossings). If one of z[i] passes
                                # zero, that is beforeEvent(z[i])*z[i] < 0, an event is triggered
                                # (allocated during instanciation according to nz).
    zDir::Vector{Int32}          # zDir[i] =  0: Root is reported for both crossing directions
                                #         =  1: Root is reported when crossing from negative to positive direction
    eweight::Vector{Float64}    # Weight vector
    last_t::Float64              # Time instant of the last call to idasol_f (residue function)
    last_norm_y::Float64         # Norm of y of the last call to idasol_f (residue function)
    hcur::MVector{1,Float64}     # Current step size argument for IDAGetCurrentStep
    order::MVector{1,Int32}       # Current order for IDAGetCurrentOrder
    r::Vector{Float64} # Residues vector used in idasol_f function
    sol_mem::SolData       # IDA pointer (to access all IDAgetXXX functions)

    y::Vector{Float64}    # Julia vector of y wrapping y_N_Vector vector
    yp::Vector{Float64}   # Julia vector of yp wrapping yp_N_Vector vector
    fulljac::Union{Matrix{Float64}, Nothing}  # Julia matrix of fulljac wrapping fulljac_N_xxx

    function IntegratorData(model::ModiaMath.AbstractSimulationModel, simulationState::ModiaMath.SimulationState)
        ny = simulationState.nx
        nz = simulationState.nz
        new(model, simulationState, simulationState.statistics,
          ones(nz), fill(0, 0), zeros(ny), typemin(Float64), 0.0,
          MVector{1,Float64}(0.0), MVector{1,Int32}(0), zeros(ny))
    end
end

include("UserData.jl")

function updateStatistics!(sol_mem, stat::ModiaMath.SimulationStatistics)
    info = 0
    info = sol_mem.steps
    stat.nSteps        += info

    info = sol_mem.jac_evals
    stat.nJac          += info

    info = sol_mem.fails
    stat.nErrTestFails += info

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
index 1 DAE as special cases). For details see [`ModiaMath.StructureOfDAE`](@ref).

During continuous integration, the DAE is solved with the
[Sundials](https://computation.llnl.gov/projects/sundials) IDA solver
(accessed via the Julia [Sundials](https://github.com/JuliaDiffEq/Sundials.jl) interface package).

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

    use_fulljac = false

    sim           = model.simulationState
    sim.model     = model
    sim.tolerance = isnan(tolerance) ? sim.defaultTolerance : convert(Float64, tolerance)
    sim.startTime = isnan(startTime) ? sim.defaultStartTime : convert(Float64, startTime)
    sim.stopTime  = isnan(stopTime)  ? sim.defaultStopTime  : convert(Float64, stopTime)
    sim.interval  = isnan(interval)  ? (isnan(sim.defaultInterval) ? (sim.stopTime - sim.startTime)/500.0 : sim.defaultInterval)  : convert(Float64, interval)
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

    if ModiaMath.isLogInfos(logger)
        println("... ModiaMath.simulate! (version ", ModiaMath.Version, " ", ModiaMath.Date, ") to simulate model: ", sim.name)
    end

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
    statistics.h0 = interval
    statistics.hMin = Inf
    statistics.hMax = 0
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
    simModel.hcur[1] = statistics.h0
    stateEvent::Bool = false
    timeEvent::Bool  = false
    flag = false
    ny = sim.nx
    nz = sim.nz

    y  = copy(init.y0)
    yp = copy(init.yp0)
    simModel.y  = y
    simModel.yp = yp
    tolAbs = tolerance * 1e-4
    if use_fulljac
        fulljac = zeros(sim.nx, sim.nx)
        simModel.fulljac = fulljac
    else
        simModel.fulljac = nothing
    end

    eventInfo = DAE.EventInfo()
    tret    = [0.0]
    jc::Int = 0  # communication point counter
    if ModiaMath.isLogInfos(logger) && nz > 0
        rootInfo = fill(0, nz)
    end

    # Create IDA
    # Allocate N_Vector storage for y and yp
    y_Vector  = deepcopy(y)
    yp_Vector = deepcopy(yp)
    relTol =  0.1 * tolerance * init.y_nominal
    for i in 1:ny
        if !init.y_errorControl[i]
            tolAbs[i] = 1e5 * init.y_nominal[i]   # switch tolerance control off
        end
    end
    # Run solution
    # Initialize solvers
    tReached = t0
    #println(tolerance, tolAbs, sim.hev, simModel.hcur[1], Inf, 6, fill(0.0, ny), 0, 0)
    mem = SolData(tolAbs, relTol, sim.hev, simModel.hcur[1], Inf, 6, fill(0.0, ny), 0, 0, 0)
    simModel.sol_mem = deepcopy(mem)
    solver = dassl()
    #println(solver)

    y0 = y_Vector
    yp0 = yp_Vector
    if use_fulljac
        # IDADlsJacFn
        jac_fun = idasol_fulljac
    end


    # Sundials.IDASetId(mem, y_states_nvector)
    #if sim.sparse
    #   nnz_jac = nnz(sim.jac)
    #   Sundials.IDAKLU(mem, ny, nnz_jac)
    #   Sundials.IDAKLUSetOrdering(mem, KLUorderingChoice)
    #   Sundials.IDASlsSetSparseJacFn(mem);
    #else
    #A  = Sundials.SUNDenseMatrix(ny, ny)
    #LS = Sundials.SUNDenseLinearSolver(y, A)
    #Sundials.IDADlsSetLinearSolver(mem, LS, A)
    #end

    # Initialize zero crossing function, if required
    hasZeroCross = nz > 0
    if hasZeroCross
        ModiaMath.DAE.getEventIndicators!(simModel.model, sim, t0, y0, yp0, simModel.z)
        simModel.zDir = zeros(nz)
        sim.zDir = simModel.zDir
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
        DAE.processEvent!(model, sim, tReached, y_Vector, yp_Vector, eventInfo)
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
            simModel.y = y
            simModel.yp = yp
            #DAE.reinitialize!(simModel.model, sim, tReached)
        end
    end

    # Loop over communication points
    cpuStartIntegration = time_ns()

    flag = false

    while tReached < tEnd && restart != ModiaMath.Terminate  #------------------ while --------------
        #println("entered")
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
        #initialize derivatives
        #flag = Sundials.IDASolve(mem, tNext, tret, simModel.y, simModel.yp, Sundials.IDA_NORMAL)
        r =  simModel.r
        tspan = (tReached, tNext)
        #println("r = ", simModel.r)

        #println("entered dassl: y = $y, yp = $yp")
        #resprob!(r, du, u, p, t) =sol_f(t, u, du, r, simModel)
        resprob!(r, du, u, p, t) = sol_f!(simModel, sim, t, u, du, r, simModel.hcur[1])
        differential_vars = [true,true,false]
        prob = DAEProblem{true}(resprob!, yp, y, tspan, differential_vars=differential_vars)
        sol = Vector{Any}
        try
            sol = solve(prob, solver, reltol = relTol*(tolAbs), abstol = tolAbs, initstep = simModel.hcur[1] )
        catch
            tolAbs *= 1000
            relTol *= 1000
            sol = solve(prob, solver, reltol = relTol, abstol = tolAbs, initstep = simModel.hcur[1] )
        end



        endInd = length(sol.t)
        simModel.yp = sol.du[end]
        if hasZeroCross

            len = length(sol.t)
            for el in 2:len
                ModiaMath.DAE.getEventIndicators!(simModel.model, sim, sol.t[el-1], sol.u[el-1], sol.du[el-1], simModel.z)
                old_z = copy(simModel.z)
                #println("old_z = $old_z")
                ModiaMath.DAE.getEventIndicators!(simModel.model, sim, sol.t[el], sol.u[el], sol.du[el], simModel.z)
                zs = old_z.*(simModel.z)
                z = copy(simModel.z)
                #println("new_z = $z, zs = $zs")
                #println("Event ind = $z")
                if (tReached!=t0)
                    flag = any(x->x<0, zs) && !flag
                    #println("flag = $flag")
                end
                if flag
                endInd = el
                println("before event at all: y =", sol.u, ", \n yp = ", sol.du)
                    break
                end
            end
        end


        nSteps = length(sol.t)-1
        if nSteps>0
            steps = zeros(nSteps)
            for it in 1:nSteps
                step = sol.t[it+1] - sol.t[it]
                steps[it] = step
            end
            simModel.hcur[1] = steps[end]*100
            statistics.h0 = steps[end]

            sim.time             = tReached
            simModel.last_t      = tReached
            simModel.last_norm_y = norm(simModel.y, Inf)
            #println("stat = ", statistics.h0)
            statistics.hMin = min(statistics.hMin, minimum(steps))
            statistics.hMax = max(statistics.hMax, maximum(steps))
        end

        yt = sol.u[endInd]
        ypt = sol.du[endInd]
        tReachedt = sol.t[endInd]
        #println("y old = $yt, yp = $ypt,  tol = $tolerance, t = $tReachedt")
        if flag
            FTOL=eps(Float64)^(1/2)
            root_found = false

            #resprob2!(r, du, u, p, t) = ModiaMath.DAE.getResidues!(simModel.model, sim, t, u, du, r, statistics.h0)
            #differential_vars = [true,true,false]
            #prob2 = DAEProblem{true}(resprob2!, yp, y, tspan, differential_vars=differential_vars)
            #sol = solve(prob2, solver, reltol = FTOL, abstol = FTOL, initstep = FTOL) #solve with better precision
            #println("sol = $sol")
            ModiaMath.DAE.getEventIndicators!(simModel.model, sim, sol.t[endInd], sol.u[endInd], sol.du[endInd], simModel.z)
            zs = old_z.*(simModel.z)
            if !any(x->x<0, zs)
                println("no root after better search, zs = $zs")
                y = sol.u[endInd]
                yp = sol.du[endInd]
                tReached = sol.t[endInd]
                flag = false
            else #binary search
                len = length(sol.u)
                ind_l = 1
                ind_r = len
                while (ind_r - ind_l > 1) #binsearch by indices
                    ind_mid = ind_l + Int32(floor((ind_r - ind_l)/2))
                    ModiaMath.DAE.getEventIndicators!(simModel.model, sim, sol.t[ind_mid], sol.u[ind_mid], sol.du[ind_mid], simModel.z)
                    if any(x->x<0, old_z.*(simModel.z))
                        ind_r = ind_mid
                    else
                        ind_l = ind_mid
                    end
                end #end binsearch by indices
                #println("indices are ind_l = $ind_l and ind_r = $ind_r")
                yMin = sol.u[ind_l]
                ypMin = sol.du[ind_l]
                tMin = sol.t[ind_l]
                yMax = sol.u[ind_r]
                ypMax = sol.du[ind_r]
                tMax = sol.t[ind_r]
                z = copy(simModel.z)
                #steps = 100
                yp = ypMax
                y = yMax
                tReached = tMax
                while !any(x->abs(x)<tolAbs^2, z) #binsearch by values approx
                #while  false #binsearch by values approx
                    tspan = (tMin, tMax)
                    #prob = DAEProblem{true}(resprob!, ypMin, yMin, tspan, differential_vars=differential_vars)
                    #sol = solve(prob, solver, reltol = relTol*tolAbs^2, abstol = tolAbs^2, initstep = simModel.hcur[1]/steps,  maxstep =  simModel.hcur[1]/steps*10 )
                    #println("sol_inside = $sol")
                    #steps = steps*10
                    yMid = (yMax + yMin)/2
                    ypMid = (ypMax + ypMin)/2
                    tMid = (tMax + tMin)/2
                    ModiaMath.DAE.getEventIndicators!(simModel.model, sim, tMid, yMid, ypMid, simModel.z)
                    z = copy(simModel.z)
                    if any(x->x<0, old_z.*(simModel.z))
                        yMax = yMid
                        ypMax = ypMid
                        tMax = tMid
                    else
                        yMin = yMid
                        ypMin = ypMid
                        tMin = tMid
                    end
                    #println("y internal = $yMid, yp = $ypMid,  z = $z, t = $tMid. amp = $amp")
                end  #end binsearch by values approx
                y = yMax
                tReached = tMax
                yp = ypMax
                #old_r = copy(simModel.r)
                #nsol_f2!(r, yp) = idasol_f(tReached, y, yp, r, simModel)
                #nsol2 = nlsolve(nsol_f2!, yp, xtol=tolAbs, ftol=tolAbs)
                #r = copy(simModel.r)
                #yp = nsol2.zero
                #println("old_r = $old_r, r = $r, zero = $nsol2")
                #if any(x->isnan(x), yp)
                #    throw(DomainError("yp is nan"))
                #end
                simModel.yp = yp

            end
        else #-- if not flag
            y = sol.u[end]
            yp = sol.du[end]
            tReached = sol.t[end]
        end #-- end if flag --

        #println("y new = $y, yp = $yp,  tol = $tolerance, t = $tReached")
        simModel.y = y
        simModel.yp = yp


        sim.time             = tReached
        simModel.last_t      = tReached
        simModel.last_norm_y = norm(simModel.y, Inf)

        stateEvent = flag
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
                println("\n y = $y, yp=$yp")
                if stateEvent
                   # Print information about the root
                    ModiaMath.DAE.getEventIndicators!(simModel.model, sim, tReached, y, yp, simModel.z)
                    rootInfo = simModel.z
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
            old_z = simModel.z
            # Trigger event
            DAE.reset!(eventInfo)
            DAE.processEvent!(model, sim, tReached, y, yp, eventInfo)
            println("\n  event = $eventInfo")
            println("\n after event y = $y, yp=$yp, t = $tReached")

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
                updateStatistics!(mem, statistics)
                simModel.y = y
                simModel.yp = yp
                DAE.reinitialize!(simModel.model, sim, tReached)
                #y = simModel.y
                #yp = simModel.yp
                ##Sundials.__IDAReInit(mem, tReached, y_N_Vector, yp_N_Vector);
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
    updateStatistics!(mem, statistics)

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
