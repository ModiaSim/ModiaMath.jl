var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#ModiaMath.jl-Documentation-1",
    "page": "Home",
    "title": "ModiaMath.jl Documentation",
    "category": "section",
    "text": "ModiaMath provides a simulation engine and other mathematical utilities for packages Modia and Modia3D that are used to model physical systems such as electrical circuits, robots, or vehicles. The recommended way is to use ModiaMath via Modia or Modia3D. However, ModiaMath is self-contained and can be also used without Modia/Modia3D.The central part of ModiaMath is a simulation engine to solve implicit index one Differential Algebraic Equations (DAEs) with and without time and state events. The theory is partially described in (Otter/Elmqvist, 2017). In particular it is shown, that a large class of DAEs can be transformed automatically to this form (including multibody systems with kinematic loops). As integrator currently IDA of the Sundials integrator suite is used (via the Julia Sundials interface package). It is planned to adapt ModiaMath to Julia package DifferentialEquations and use IDA and other appropriate integrators via this package in the future.Additionally, ModiaMath provides functions to perform plotting in a convenient way, to generate and use rotation matrices and quaternions for kinematic transformations in 3D, and to provide an infrastructure for DAE variables as needed by Modia3D."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "The package is registered in METADATA.jl and can be installed in the following way (Julia >= 1.0 is required):julia> ]add ModiaMathModiaMath uses PyPlot for plotting. If PyPlot is not available in your current Julia environment an information message is printed and all ModiaMath.plot(..) calls are ignored. In order that plot windows are displayed, you need to add PyPlot to your current environment via Pkg.add(\"PyPlot\"). Often this automatic installation fails and it is recommended to follow instead the instructions Installing PyPlot in a robust way."
},

{
    "location": "index.html#Release-Notes-1",
    "page": "Home",
    "title": "Release Notes",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Version-0.4.0-1",
    "page": "Home",
    "title": "Version 0.4.0",
    "category": "section",
    "text": "The first version that requires Julia >= 1.0 (all Julia 0.6 and 0.7 code was removed).ModiaMath.SimulationState has a new keyword argument structureOfDAE. It is now possible to define that the DAE is linear in all derivatives or all derivatives are explicitely solved (= ODE). If one of the new structures is selected,  initialization and re-initialization is more robust and more efficient (in particular no implicit Euler step is used).ModiaMath.plot supports variables that can be computed from a vector of struct result (for details, see docu of plot).Using newest versions of all used packages (since issues with Sundials).Docu of SimulationState and of StructureOfDAE introduced.Version information adapted to Modia style (Version/Date constants).Info message of simulation start only printed if log=true."
},

{
    "location": "index.html#Version-0.3.1-1",
    "page": "Home",
    "title": "Version 0.3.1",
    "category": "section",
    "text": "Initialization issue corrected\nUse PyPlot.pygui(true) when importing ModiaMath in order that separate plots windows are used and no inline plots, independent of the used environment (especially made so that plots in vs-code are not inlined in the vs-code GUI).\nModiaMath.plot supports now also Vectors of structs. E.g. Vector{ThermodynamicState} and plot(result, \"state.p\") is then possible, where state is an instance of the struct and \"p\" is a fieldname with a scalar Number value."
},

{
    "location": "index.html#Version-0.3.0-1",
    "page": "Home",
    "title": "Version 0.3.0",
    "category": "section",
    "text": "ModiaMath initialization and re-initialization improved:\nInitialization and re-initialization considerably changed (now only one nonlinear equation is solved and no longer two nonlinear equations).\nBy this change, several issues in the initialization were fixed.\nx_nominal also supported for ModiaMath.Variables."
},

{
    "location": "index.html#Version-0.2.6-1",
    "page": "Home",
    "title": "Version 0.2.6",
    "category": "section",
    "text": "ModiaMath initialization and simulation improved:\nIssue with defaultInterval fixed: Previously: defaultInterval = (defaultStopTime - defaultStartTime)/500.0. However, if stopTime was set to 10e5, then the defaultInterval was much too small. This was changed so that defaultInterval is either explicitly given, or it is computed from the actual values of StopTime and StartTime.\nx_nominal introduced in constructor of simulation state (DAE.SimulationState) and of Modia interface (ModiaToModiaMath.ModiaSimulationModel): If x_nominal is explicitly provided in these constructor calls, it is used. Otherwise x_nominal[i] is set to max(abs(x_start[i]), 1e-7).\nx_nominal is newly used in KINSOL to scale the unknowns and is used to compute the absolute tolerance for Sundials IDA.\nA bug was corrected in the initialization function, where rScale was reported to KINSOL, but was also used for scaling of the residue (although this is performed in KINSOL).\nIf log=true, name-, start-, fixed-, nominal-values of the x-vector are printed before initialization starts.New function ModiaMath.solveOneNonlinearEquation: Determines the solution of one non-linear algebraic equation y=f(u) in one unknown u in a reliable and efficient way (using Brents algorithm)."
},

{
    "location": "index.html#Version-0.2.5-1",
    "page": "Home",
    "title": "Version 0.2.5",
    "category": "section",
    "text": "ModiaMath initialization improved:\nKINSOL tolerance FTOL takes IDA tolerance into acount.\nMore information about the model is provided in case KINSOL fails.Interface to Modia improved: New keyword option hev (step size used for implicit Euler during initialization and at events) added to \"simulate\".ModiaMath result handling improved:\nThe result may contain single values (such as parameter J=0.1). When using plot(..) on such a result variable, a constant line is plotted between the first and the last point of the x-axis.\nNew function ModiaMath.resultTable(result) to transform the variable information in the result data structure in a DataFrames table (containing variable name + type + size ) that can then be printed.\nNew functions ModiaMath.closeFigure(figure) and ModiaMath.closeAllFigures() to close a specific figure or close all figures.\nChanging some default options of PyPlot. In particular, (a) the labels on the x- and y-axis use exponential notation (e.g. 1e5), if the numbers are larger as 1e4 or smaller as 1e-3 (PyPlot default is 1e7 and 1e-7), (b) smaller fonts and linewidth are used. Default options are changed via PyCall. If PyCall is not in the Julia environment, PyCall is added.\nIf a vector or matrix of subplots is defined, then the x-axis labels are only displayed for the last subplot row.\nResult data is always converted to Float64 before passing it to PyPlot. Therefore, Julia numbers can be plotted, even if the types of the numbers are not supported by PyPlot (for example rational numbers can be plotted)."
},

{
    "location": "index.html#Version-0.2.4-1",
    "page": "Home",
    "title": "Version 0.2.4",
    "category": "section",
    "text": "Non-backwards compatible change of ModiaMath.plot(..): The function was changed to only support result dictionaries of the type Dict{AbstractString,Any} to simplify implementation and maintenance and to use the identical dictionary type as in Modia.ModiaMath.plot(..) improved:\nFor the dictionary type used by Modia, all elements of a variable vector are plotted by giving the variable name, or one specific element of the variable vector is plotted by giving the name and the vector index.\nA new keyword argument maxLegend=10 introduced: If the number of entries in a legend exceeds this number, no legend is included in the plot. Note, the curves can still be identified by clicking in the tool bar of the plot window on button Edit axis, curve ...\nIf the plot function is used in the wrong way (e.g. signals shall be plotted that are not in the result dictionary, or selecting a variable vector as x-axis), warning messages are printed and the call is ignored.Dependent packages updated to their newest versions. Especially, warnings from Unitful do no longer occur, due to the update to version 0.12.0.Issue with PyPlot warnings fixed: When ModiaMath was used in another package (e.g. in Modia or Modia3D), strange warning messages appeared when importing it. This issue was fixed with the solution sketched in discourse."
},

{
    "location": "index.html#Version-0.2.3-1",
    "page": "Home",
    "title": "Version 0.2.3",
    "category": "section",
    "text": "Wrong UUID of ModiaMath corrected (did not correspond to the UUID in Julias METADATA)."
},

{
    "location": "index.html#Version-0.2.2-1",
    "page": "Home",
    "title": "Version 0.2.2",
    "category": "section",
    "text": "PyPlot was removed from the REQUIRE and Project.toml files and code was added, so that PyPlot is automatically imported in ModiaMath if it is available in the current environment of the user. The benefit is that ModiaMath can be used, even if PyPlot is not installed. This is especially useful for ContinuousIntegration, because automatic installation of PyPlot often fails. Inspect the wiki page Installing PyPlot in a robust way to install PyPlot in a robust way.All extra packages used in examples and tests are now referenced via ModiaMath (for example using ModiaMath.StaticArrays instead of using StaticArrays). The benefit is that all examples and tests can be directly executed with include (for example: import ModiaMath; include($(ModiaMath.path)/examples/Simulate_Pendulum.jl)) provided ModiaMath is in the current environment. Previously, it was assumed that these extra packages are present in the users environment and an error occured, if this was not the case.New arguments prefix and reuse added to ModiaMath.plot(..) to add new plots (e.g. from a new simulation run) to an existing figure, without clearing the figure beforehand."
},

{
    "location": "index.html#Version-0.2.1-1",
    "page": "Home",
    "title": "Version 0.2.1",
    "category": "section",
    "text": "Adapted to Julia 0.7 and 1.0 (including using new package manager via Project.toml, Manifest.toml files and adapting the README.md and documentation files).Buttons for Travis CL, code coverage and docs added to README.md file. (this includes adaptations to .travis.yml file). Note, Appveyor CL is not yet activated (although appveyor.yml file is present).Interface call to KINSOL slightly improved to avoid a gc-crash.In case KINSOL fails, more run-time information is added to the error message.Indentation was changed consistently to 4 spaces."
},

{
    "location": "index.html#Version-0.2.0-1",
    "page": "Home",
    "title": "Version 0.2.0",
    "category": "section",
    "text": "First public release (for Julia 0.6)"
},

{
    "location": "man/Overview.html#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "man/Overview.html#Overview-1",
    "page": "Overview",
    "title": "Overview",
    "category": "section",
    "text": ""
},

{
    "location": "man/Overview.html#Mathematical-Description-1",
    "page": "Overview",
    "title": "Mathematical Description",
    "category": "section",
    "text": "ModiaMath provides a simulation engine and other mathematical utilities  to solve initial value problems of the following implicit index 1 DAE with events (containing an ODE or a semi-explicit index 1 DAE as special cases; x(t), z(t) are real vectors):beginalign\r\n z = f_z(x t) \r\n 0 = f_d(dotx x t z_i  0) \r\n 0 = f_c(x t z_i  0) \r\n J = left fracpartial f_dpartial dotx  \r\n             fracpartial f_cpartial x right  textis regular\r\nendalignwherelim_epsilon rightarrow 0 x(t_0 - epsilon) = x_0^-is given. Equations z=z(t) are zero-crossing functions. Whenever a z_i(t) crosses zero  the integration is halted, functions f_d f_c might be changed and  afterwards integration is restarted. If the Jacobian J is regular, the DAE has an index 1 (= by differentiating f_c once, the system can be transformed to an ODE).Given a DAE model, ModiaMath assumes that J is regular for all time instants. If this condition is violated, initialization and simulation will usually fail and an error message of the form \"Solver does not converge\" might appear. Note, ModiaMath does not check this condition and can therefore not provide better diagnostics in such cases.Initial conditions x_0^- must be provided before simulation can start. They need  not to fulfil the constraint equations, so f_c (x_0^-t_0 )0 is allowed. If this is the case, initialization will simulate for an infinitesimal small time instant  so that x_0^- changes discontinuously to x_0^+ with f_c (x_0^+t_0 )=0.  Note, dotx is a Dirac impulse if x changes discontinuously at initialization.As shown in xxx, every DAE can be transformed to the form above, at least in principal. In yyy algorithms are proposed to automatically transform a large class of DAEs to this form without solving algebraic equations and retaining the sparsity of the equations.  This may require to analytically differentiating equations. The algorithms of this paper are implemented in the Julia package Modia which in turn uses ModiaMath. In Modia3D the transformation to this form is built into the package itself.Note, the above DAE could be further transformed to an ODE (dotx = f(xt)), but then the evaluation of function f(xt) might require to solve local linear and/or nonlinear equation systems. Furthermore, there are systems  (for example ModiaMath/examples/Simulate_FreeBodyRotation.jl) where the ODE states x need to be dynamically changed during simulation.It is highly recommended to use Modia or Modia3D for simulating DAEs because this is much simpler and less error prone as when utilizing ModiaMath directly. However, ModiaMath can be also used without Modia or Modia3D. In this case, basically one Julia function with the following interface has to be provided (and in this function specific utility functions can be called)getModelResidues(m::AbstractSimulationModel, t::Float64, x::Vector{Float64},  \r\n                 derx::Vector{Float64}, r::Vector{Float64}where r is the vector of residues (r = left f_d f_c right). Given the  simulation model m (= a mutable struct), the actual time instant t, the DAE variables x(t) and their derivatives derx(t), the function has to compute the residue vector r(t). In directory ModiaMath/examples/withoutMacros_withoutVariables several examples are present that are based on this interface.In order to simplify the definition of direct ModiaMath models (to evaluate and test ModiaMath  functionality), the macro @component has been introduced. The examples in directory ModiaMath/examples/xxx.jl use this model definition. In directory ModiaMath/examples/withoutMacros/xxx.jl the same examples are present, however, the macro has been manually expanded (to show and test the result of the macro). The @component does not yet support events. If events are present in a model, the model has to be defined as shown in the examples of directory ModiaMath/examples/withoutMacros_withoutVariables."
},

{
    "location": "man/Overview.html#Getting-Started-1",
    "page": "Overview",
    "title": "Getting Started",
    "category": "section",
    "text": "You can just past the following code into the Julia REPL."
},

{
    "location": "man/Overview.html#To-define-a-model-1",
    "page": "Overview",
    "title": "To define a model",
    "category": "section",
    "text": "(note, it is simpler and less error prone to define a model with Modia or Modia3D):  using ModiaMath\r\n  using StaticArrays\r\n\r\n  @component Pendulum(;L=1.0, m=1.0, d=0.1, g=9.81) begin\r\n     phi = RealScalar(start=pi/2, unit=\"rad\"    , fixed=true,               numericType=ModiaMath.XD_EXP)\r\n     w   = RealScalar(start=0.0 , unit=\"rad/s\"  , fixed=true, integral=phi, numericType=ModiaMath.XD_EXP)\r\n     a   = RealScalar(            unit=\"rad/s^2\",             integral=w  , numericType=ModiaMath.DER_XD_EXP) \r\n     r   = RealSVector{2}(        unit=\"m\"      ,                           numericType=ModiaMath.WC)\r\n  end;\r\n\r\n  function ModiaMath.computeVariables!(p::Pendulum, sim::ModiaMath.SimulationState)  \r\n     L = p.L; m = p.m; d = p.d; g = p.g; phi = p.phi.value; w = p.w.value\r\n   \r\n     p.a.value = (-m*g*L*sin(phi) - d*w) / (m*L^2)\r\n\r\n     if ModiaMath.isStoreResult(sim)\r\n        p.r.value = @SVector [L*sin(phi), -L*cos(phi)]\r\n     end\r\n  end;\r\n\r\n  simulationModel = ModiaMath.SimulationModel(Pendulum(L=0.8, m=0.5, d=0.2), stopTime=5.0);\r\n"
},

{
    "location": "man/Overview.html#To-simulate-a-model-and-plot-results:-1",
    "page": "Overview",
    "title": "To simulate a model and plot results:",
    "category": "section",
    "text": "  result = ModiaMath.simulate!(simulationModel; log=true);\r\n  plot(result, [(:phi, :w) :a])This results in:<img src=\"../../resources/images/pendulumPlot.svg\">"
},

{
    "location": "man/Overview.html#To-run-examples-and-tests:-1",
    "page": "Overview",
    "title": "To run examples and tests:",
    "category": "section",
    "text": "# run examples\r\nimport ModiaMath\r\ninclude(\"$(ModiaMath.path)/examples/Simulate_Pendulum.jl\")         # ODE as index-0 DAE\r\ninclude(\"$(ModiaMath.path)/examples/Simulate_FreeBodyRotation.jl\") # index-1 DAE\r\ninclude(\"$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_PendulumDAE.jl\") # index-3 DAE\r\ninclude(\"$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_SimpleStateEvents.jl\")\r\ninclude(\"$(ModiaMath.path)/examples/withoutMacros_withoutVariables/Simulate_BouncingBall.jl\")\r\n\r\n# run all tests\r\ninclude(\"$(ModiaMath.path)/test/runtests.jl\")"
},

{
    "location": "man/Overview.html#Package-structure-1",
    "page": "Overview",
    "title": "Package structure",
    "category": "section",
    "text": "The ModiaMath package is structured in the a set of sub-modules. The most important are:ModiaMath.SimulationEngine\nThe engine to simulate implicit index 1 DAEs with events.ModiaMath.DAE\nInterface between the SimulationEngine and the index 1 DAE model (e.g. initialization and event iteration is performed here).ModiaMath.Result\nThe plot function of this module allows to plot the result data of the simulation engine by giving the signal names. With tuples and/or vectors/matrices of signal names, the window layout of the figures is defined. The legends/labels of the plots are automatically constructed by the signal names and their units.ModiaMath.Variables\nProvides Variable types to define properties of the variables on a higher level and copy automatically the interface vectors from the integrator into the variables and vice versa.ModiaMath.Frames\nFunctions that generate and operate on frames, that is coordinate systems in 3D. The orientation of a frame is described either with a 3x3 rotation matrix or with a  quaternion vector. This module is currently mainly used from Modia3D, but the functionality is useful for all 3D programs."
},

{
    "location": "man/Plans.html#",
    "page": "Plans for version 1.0",
    "title": "Plans for version 1.0",
    "category": "page",
    "text": ""
},

{
    "location": "man/Plans.html#Plans-for-version-1.0-1",
    "page": "Plans for version 1.0",
    "title": "Plans for version 1.0",
    "category": "section",
    "text": "ModiaMath is not yet ready and should not be used for production simulation runs. The following features are planned to be implemented before ModiaMath version 1.0:General\nSupport of all features needed by Modia 1.0 and Modia3D 1.0.\nImproved documentation.SimulationEngine\nSparse Jacobian for IDA (mainly available, but not yet transferred to restructured simulation engine).\nSparse Jacobian for nonlinear solver during initialization.\nUser-supplied dense Jacobian for IDA (pre-requisite for next step).\nSupport for mue-stabilization technique for DAEs with index > 1 (using Jacobian information).\nChange basic engine from IDA to package DifferentialEquations (and access IDA via this interface).\nIf model is an ODE, support ODE integrators from DifferentialEquations package.\nNew inquiry function isCompletedIntegratorStep(..) to inquire when the integrator step  is completed. This function shall have the same functionality as fmi2CompletedIntegratorStep and requires to restructure the integrator while-loop since the integrator must return after every completed step.\nLinearization of the DAE after initialization and at the actual integration time instant in order that linear analysis and synthesis methods can be applied.NonlinearEquations\nTransfering the FORTRAN code hybrd.f from Minpack. to a native Julia implementation. It is expected that this solver is more robust as Sundials KINSOL (implemented with adaptations, but not yet fully tested).\nUsing Minpack during initialization.Result\nInteger and Boolean result variables.\nSeveral time vectors for different result data.\nParameters (time vector has length 2).\nClocked variables (time vector has only values at clock ticks).\nAutomatic plot layout for clocked variables.\nSupport of another plot package (besides PyPlot). Probably best to  support Plots.\nUse named tuples to define line style elements, e.g.\nplot(result, ((name=:phi1, color=:blue), (name=:phi2, color=:red)) )Frames\nSupport splines for the interpolation of Frames (currently only linear interpolation supported)."
},

{
    "location": "lib/SimulationEngine.html#",
    "page": "SimulationEngine",
    "title": "SimulationEngine",
    "category": "page",
    "text": ""
},

{
    "location": "lib/SimulationEngine.html#ModiaMath.SimulationEngine",
    "page": "SimulationEngine",
    "title": "ModiaMath.SimulationEngine",
    "category": "module",
    "text": "module ModiaMath.SimulationEngine\n\nSimulation engine for implicit index 1 DAE models with events.\n\nMain developer\n\nMartin Otter,  DLR - Institute of System Dynamics and Control\n\n\n\n\n\n"
},

{
    "location": "lib/SimulationEngine.html#ModiaMath.SimulationEngine.simulate!-Tuple{ModiaMath.AbstractSimulationModel}",
    "page": "SimulationEngine",
    "title": "ModiaMath.SimulationEngine.simulate!",
    "category": "method",
    "text": "ModiaMath.simulate!(simulationModel; log=false, startTime=NaN, stopTime=NaN, \n                                                tolerance=NaN, interval=NaN)\n\nSimulates a DAE simulationModel that is defined with package Modia, package Modia3D or with the ModiaMath.@component macro. The DAE is mathematically described as  implicit, index 1 DAE with events (containing an ODE or a semi-explicit  index 1 DAE as special cases):\n\nbeginalign\n 0       = f_d(dotx x t z_pos) \n 0       = f_c(x t z_pos) \n z       = f_z(xt) \n z_pos = textevent()    z  0    textlast(z_pos) \n J       = left fracpartial f_dpartial dotx  \n                    fracpartial f_cpartial x right  textis regular\nendalign\n\nwith initial conditions x_0^-:\n\nlim_epsilon rightarrow 0 x(t_0 - epsilon) = x_0^-\n\nDuring continuous integration, equation system (1)-(4) is solved with the Sundials IDA solver (accessed via the Julia Sundials interface package). ModiaMath assumes that J (5) is regular for all time instants. If this condition is violated, initialization and simulation will usually fail and an error message of the form \"Solver does not converge\" might appear. Note, ModiaMath does not check this condition and can therefore not provide better diagnostics in such cases.\n\nIf one of the elements of z crosses zero, an event is triggered and simulation is halted. At an event, equation z_pos = z  0 (element wise) is added. The equation system (1)-(4) is then solved with a fixed-point iteration scheme (= event iteration). Afterwards, integration is  restarted and z_pos keeps its value until the next event occurs.\n\nInitial conditions x_0^- must be provided before simulation can start.  Best is if they fulfil the constraint equation 0 = f_c(x_0^- t_0 z  0). If this is not the case, initialization will simulate for an infinitesimal small time instant  so that x_0^- changes discontinuously to x_0^+ with f_c (x_0^+ t_0 z  0 )=0.  Note, dotx is a Dirac impulse in this case.\n\nInput arguments\n\nsimulationModel::ModiaMath.AbstractSimulationModel: Model struct (generated with Modia, Modia3D or ModiaMath.@component).\nlog::Bool: = true, if logging is enabled, otherwise it is disabled.\nstartTime::Float64: Start time of the simulation in [s].                       If startTime=NaN, the default startTime is used that is defined by the simulationModel.\nstopTime::Float64: Stop time of the simulation in [s].                       If stopTime=NaN, the default stopTime is used that is defined by the simulationModel.\ntolerance::Float64: The relative tolerance for the integration. The absolute tolerance is computed as 0.1*tolerance*nominal(variable) where nominal(variable) is the nominal value of the variable.  If tolerance=NaN, the default tolerance is used that is defined by the simulationModel.\ninterval::Float64: Output interval for results in [s]. If events occur, the event time instants are additionally added to the result. If interval=NaN, the default interval is used that is defined by the simulationModel.\n\n\n\n\n\n"
},

{
    "location": "lib/SimulationEngine.html#SimulationEngine-1",
    "page": "SimulationEngine",
    "title": "SimulationEngine",
    "category": "section",
    "text": "Modules = [ModiaMath.SimulationEngine]\r\nPrivate = false\r\nOrder   = [:module, :type, :macro, :function]"
},

{
    "location": "lib/Result.html#",
    "page": "Result",
    "title": "Result",
    "category": "page",
    "text": ""
},

{
    "location": "lib/Result.html#ModiaMath.Result",
    "page": "Result",
    "title": "ModiaMath.Result",
    "category": "module",
    "text": "module ModiaMath.Result\n\nOrganize and plot simulation result data (time series). The result data of ModiaMath.simulate! is returned in one of the formats supported by this module. The ModiaMath.plot function of this module allows to plot the result data by giving the signal names. The legends/labels of the plots are automatically constructed by the signal names and their unit. Example:\n\nModiaMath.plot(result, [ (:phi,:r)      (:phi,:phi2,:w);\n                         (:w,:w2,:phi2) (:phi,:w)      ], \n               heading=\"Matrix of plots\")\n\ngenerates the following plot:\n\n(Image: Matrix-of-Plots)\n\nMain developer\n\nMartin Otter,  DLR - Institute of System Dynamics and Control\n\n\n\n\n\n"
},

{
    "location": "lib/Result.html#ModiaMath.Result.RawResult",
    "page": "Result",
    "title": "ModiaMath.Result.RawResult",
    "category": "type",
    "text": "mutable struct RawResult\n\nHold result data in a Float64 matrix (first column time, i-th column variable i). Typical usage:\n\n    raw = RawResult(100,3) # Initial storage: 100 time points, 3 variables\n    storeRawResult!(raw, [1.0, 2.0, 1.0])\n    storeRawResult!(raw, [2.0, 3.0, 5.0])\n      ...\n    res = getDictResult(res, [\"time\", \"r[1]\", \"r[2]\"]\n    \n    plot(res[\"time\"], res[\"r[1]\"])\n\nIf initial storage is not sufficient, it is automatically doubled.\n\n\n\n\n\n"
},

{
    "location": "lib/Result.html#ModiaMath.Result.ResultWithVariables",
    "page": "Result",
    "title": "ModiaMath.Result.ResultWithVariables",
    "category": "type",
    "text": "mutable struct ResultWithVariables\n\nStruct that is generated as a result of ModiaMath.simulate!(...), if the model is a Modia3D.AbstractComponentWithVariables struct.\n\n\n\n\n\n"
},

{
    "location": "lib/Result.html#ModiaMath.Result.closeAllFigures-Tuple{}",
    "page": "Result",
    "title": "ModiaMath.Result.closeAllFigures",
    "category": "method",
    "text": "ModiaMath.closeAllFigures()\n\nCloses all open figures.\n\n\n\n\n\n"
},

{
    "location": "lib/Result.html#ModiaMath.Result.closeFigure-Tuple{Int64}",
    "page": "Result",
    "title": "ModiaMath.Result.closeFigure",
    "category": "method",
    "text": "ModiaMath.closeFigure(figure::Int)\n\nCloses figure with figure number figure.\n\n\n\n\n\n"
},

{
    "location": "lib/Result.html#ModiaMath.Result.getStringDictResult-Tuple{Any,ModiaMath.Result.RawResult}",
    "page": "Result",
    "title": "ModiaMath.Result.getStringDictResult",
    "category": "method",
    "text": "getStringDictResult(model, res)\n\nReturn dictionary of result, from raw result data structure raw.\n\n\n\n\n\n"
},

{
    "location": "lib/Result.html#ModiaMath.Result.plot-Tuple{Any,Symbol}",
    "page": "Result",
    "title": "ModiaMath.Result.plot",
    "category": "method",
    "text": "ModiaMath.plot(result, names; heading=\"\", grid=true, xAxis= :time, \n               figure=1, prefix=\"\", reuse=false, maxLegend=10)\n\nPlot time series of the result defined by the names keys (Symbol or String). The keys (and their units, if available in the result) are automatically used as legend. Units can be either added by using package Unitful if result is just a dictionary, or it can be added by using package ModiaMath.Result, where units are defined as elements of the variable definition. \n\nArguments\n\nArgument result maybe one of the following:\n\nA dictionary Dict{AbstractString,Any}. The Dict Value can be \n(a) Vector{Number},\n(b) Matrix{Number},\n(c) a vector of structs where the field that shall be plotted is a Number, or\n(d) a vector of structs and a function is provided to compute a Number from the struct (this function must be stored as value in a dictionary Dict{AbstractString,Function} and this dictionary is returned from function ModiaMath.variablesDependingOnStruct(struct); see last example below).\nNote, before passing data to the plot package, it is converted to Float64. This allows to, for example, also plot rational numbers, even if not supported by the plot package. \nAn instance of struct ModiaMath.Result\nAn object for which function ModiaMath.resultTimeSeries is defined.\n\nArgument names defines the diagrams to be drawn and the time series to be included in the respective diagram: \n\nIf names is a Symbol or String, generate one diagram with one time series.\nIf names is a Tuple of Symbols/Strings, generate one diagram with the time series of the given keys\nIf names is a Vector or Matrix of Symbols/Strings/Tuples, generate a vector or matrix of diagrams.\n\nRemaining arguments:\n\nheading::AbstractString: Optional heading above the diagram.\ngrid::Bool: Optional grid.\nxAxis: Name of x-axis (Symbol or AbstractString).\nfigure::Int: Integer identifier of the window in which the diagrams shall be drawn.\nprefix::AbstractString: String that is appended in front of every legend label (useful especially if reuse=true)\nreuse::Bool: If figure already exists and reuse=false, clear the figure before adding the plot.\nmaxLegend::Int: If the number of legend entries in one plot command > maxLegend, the legend is suppressed. All curves have still their names as labels. The curves can be inspected by their names by clicking in the toolbar of the plot on button Edit axis, curve .. and then on Curves.\n\nExamples\n\nimport ModiaMath\nusing Unitful\n\nt = linspace(0.0, 10.0, 100)\nresult = Dict{AbstractString,Any}(\n            \"time\" => t*u\"s\", \"phi1\" => sin.(t)u\"rad\"  , \"phi2\" => 0.5*sin.(t),\n                              \"w1\"   => cos.(t)u\"rad/s\", \"w2\"   => 0.6*cos.(t))\n\n# 1 signal in one diagram\n#   (legend = \"phi1 [rad]\")\nModiaMath.plot(result, :phi1)   \n\n# 3 signals in one diagram                                 \nModiaMath.plot(result, (\"phi1\", :phi2, :w1), figure=2)\n\n# 3 diagrams in form of a vector (every diagram has one signal)                 \nModiaMath.plot(result, [:phi1, :phi2, :w1], figure=3)     \n\n# 4 diagrams in form of a matrix (every diagram has one signal)          \nModiaMath.plot(result, [\"phi1\" \"phi2\";\n                        \"w1\"   \"w2\"   ], figure=4)     \n\n# 2 diagrams in form of a vector           \nModiaMath.plot(result, [ (:phi1,:phi2), (:w1) ], figure=5)           \n\n# 4 diagrams in form of a matrix\nModiaMath.plot(result, [ (:phi1,)           (:phi2,:w1);\n                         (:phi1,:phi2,:w1)  (:w2,)     ],figure=6)  \n\n# Plot w1=f(phi1) in one diagram \nModiaMath.plot(result, :w1, xAxis=:phi1, figure=7)    \n\n# Append signal of the next simulation run to figure=1\n# (legend = \"Sim 2: phi1 [rad]\")\nresult[:phi1] = 0.5*result[:phi1]\nModiaMath.plot(result, :phi1, prefix=\"Sim 2: \", reuse=true)\n\n# Compute and plot variables that depend on the result vector\nmutable struct MyThermodynamicState\n    p::Float64\n    T::Float64\nend\nspecificEnthalpy(state::MyThermodynamicState) = 2.0*state.T\ndynamicViscosity(state::MyThermodynamicState) = 2.0*state.p\n\nconst dependentVariables = Dict{AbstractString,Function}(\"h\"   => specificEnthalpy,\n                                                         \"eta\" => dynamicViscosity)\n\nModiaMath.variablesDependingOnStruct(state::MyThermodynamicState) = dependentVariables\n\nstate = MyThermodynamicState[]\nfor i = 1:length(t)\n    push!(state, MyThermodynamicState(i*2.0,i*3.0))\nend\nresult = Dict{AbstractString,Any}(\"time\" => t, \"state\" => state)\n\nModiaMath.plot(result, (\"state.p\", \"state.T\", \"state.h\", \"state.eta\"))\n\nThe 5th example above (2 diagrams in form of a vector) give the following plot:\n\n(Image: Figure 5)\n\n\n\n\n\n"
},

{
    "location": "lib/Result.html#ModiaMath.Result.resultTable-Tuple{Dict{AbstractString,Any}}",
    "page": "Result",
    "title": "ModiaMath.Result.resultTable",
    "category": "method",
    "text": "table = resultTable(result)\n\nReturn the variables stored in result in form of a DataFrames table (which can then be printed/showed in various forms).\n\nBase.show(io, result) is defined to print resultTable(result), in case result is of type ModiaMath.ResultWithVariables.\n\nExamples\n\nimport ModiaMath\nusing  Unitful\n\nt = range(0.0, stop=10.0, length=100)\nresult = Dict{AbstractString,Any}()\nresult[\"time\"] = t * u\"s\";\nresult[\"phi\"]  = sin.(t)u\"rad\";\n\n# Print table of the variables that are stored in result\nprintln(\"result variables = \", ModiaMath.resultTable(result))\n\n# Results in\nresult variables =\n│ Row │ name   │ elType  │ sizeOrValue   │ unit   │\n│     │ String │ String  │ String        │ String │\n├─────┼────────┼─────────┼───────────────┼────────┤\n│ 1   │ phi    │ Float64 │ (100,)        │ rad    │\n│ 2   │ time   │ Float64 │ (100,)        │ s      │\n\n\n\n\n\n\n"
},

{
    "location": "lib/Result.html#ModiaMath.Result.resultTimeSeries-Tuple{Dict{AbstractString,Any},Any,Bool,Any}",
    "page": "Result",
    "title": "ModiaMath.Result.resultTimeSeries",
    "category": "method",
    "text": "(xsig, xsigLegend, ysig, ysigLegend) = \n       ModiaMath.resultTimeSeries(result, name, xLabel::Bool, xAxis)\n\nFor a desired result data structure, this function has to be provided to return the x-vector (xsig), the y-vector (ysig) and the legend of the x-vector (xsigLegend), and of the y-vector (ysigLegend) as Strings, given  the key of the y-vector (name) and the key of the x-vector (xAxis). If xLabel=false the legend of the x-vector should be an empty string (\"\").\n\n\n\n\n\n"
},

{
    "location": "lib/Result.html#ModiaMath.Result.storeRawResult!-Tuple{ModiaMath.Result.RawResult,Array{Float64,1}}",
    "page": "Result",
    "title": "ModiaMath.Result.storeRawResult!",
    "category": "method",
    "text": "storeRawResult!(res, v::Vector{Float64})\n\nStore vector v in result data structure res.\n\n\n\n\n\n"
},

{
    "location": "lib/Result.html#Result-1",
    "page": "Result",
    "title": "Result",
    "category": "section",
    "text": "Modules = [ModiaMath.Result]\r\nPrivate = false"
},

{
    "location": "lib/DAE.html#",
    "page": "DAE",
    "title": "DAE",
    "category": "page",
    "text": ""
},

{
    "location": "lib/DAE.html#ModiaMath.DAE",
    "page": "DAE",
    "title": "ModiaMath.DAE",
    "category": "module",
    "text": "module ModiaMath.DAE\n\nInterface between the ModiaMath.SimulationEngine and the index 1 DAE model. A DAE model is a struct that has a required field  simulationState::ModiaMath.SimulationState in which the main properties of the DAE model are reported to the simulation engine:\n\n# DAE model ModelName\nmutable struct ModelName <: ModiaMath.AbstractSimulationModel\n    simulationState::ModiaMath.SimulationState\n\n    # other definitions (e.g. parameters of model)\nend\n\nThe following functions can be called in the DAE model to inquire information about the simulation state:\n\nModiaMath.getTime\nModiaMath.getStartTime\nModiaMath.getStopTime\nModiaMath.getTolerance\nModiaMath.isInitial\nModiaMath.isTerminal\nModiaMath.isEvent\nModiaMath.isZeroCrossing\nModiaMath.isAfterSimulationStart\nModiaMath.isStoreResult\nModiaMath.isLogInfos\nModiaMath.isLogWarnings\nModiaMath.isLogEvents\n\nThe following functions can be called in the DAE model to set properties in the simulation engine:\n\nModiaMath.setNominal!\nModiaMath.setNextEvent!\nModiaMath.positive!\nModiaMath.negative!\nModiaMath.change!\nModiaMath.edge!\n\nThe following functions can be either called in the DAE model or they can be called on a simulation model (before or after ModiaMath.simulate!(simulationModel, ...) is called).\n\nModiaMath.logOn!\nModiaMath.logOff!\nModiaMath.setLogCategories!\n\nMain developer\n\nMartin Otter,  DLR - Institute of System Dynamics and Control\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.EventRestart",
    "page": "DAE",
    "title": "ModiaMath.DAE.EventRestart",
    "category": "type",
    "text": "@enum EventRestart NoRestart Restart FullRestart Terminate\n\nDefine how to continue or restart integration after an event. Usually, Restart should be used. Only in special cases, the other flags are useful.\n\nNoRestart, continue integration without restarting the integrator\nRestart, restart integrator\nFullRestart, restart integrator and optionally exchange simulationState (so dimensions may change)\nTerminate, terminate integration\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.SimulationState",
    "page": "DAE",
    "title": "ModiaMath.DAE.SimulationState",
    "category": "type",
    "text": "simulationState = SimulationState(name, getModelResidues!, x_start, \n                                  getVariableName; kwargs...)\n\nReturn a simulationState object that is described by a DAE with one of the supported structures defined with enumeration StructureOfDAE.\n\nA model that shall be simulated with function ModiaMath.simulate!(model, ...) is required to be defined as:\n\nmutable struct ModelName <: ModiaMath.AbstractSimulationModel\n    simulationState::ModiaMath.SimulationState\n\n    # other definitions (e.g. parameters of model)\nend\n\nKeyword arguments defaults\nstructureOfDAE ImplicitIndexOneDAE (see StructureOfDAE)\nnc 0\nnz 0\nnw 0\nzDir fill(0, nz)\nx_fixed fill(false, length(x_start))\nx_nominal fill(NaN, length(x_start))\nx_errorControl fill(true, length(x_start))\njac nothing\nmaxSparsity 0.1\nhev 1e-8\nscaleConstraintsAtEvents true,\ngetResultNames ModiaMath.getResultNames\nstoreResult! ModiaMath.storeRawResult!\ngetResult ModiaMath.getStringDictResult\ndefaultTolerance 1e-4\ndefaultStartTime 0.0\ndefaultStopTime 1.0\ndefaultInterval NaN\n\nRequired arguments\n\nname::Union{AbstractString,Symbol}: Name of model\ngetModelResidues!::Function: Function with arguments (model,t,x,derx,r,w) to compute the residues r and auxiliary variables w from time t, vector x and its time derivative derx. \n\nx_start::Vector{Float64}: Start values of x.\ngetVariableName::Function=ModiaMath.defaultVariableName: Function that returns the name of a variable, given its type and its index.\n\nOptional (keyword) arguments:\n\nstructureOfDAE::StructureOfDAE: Structure of DAE. Try to not use the default ImplicitIndexOneDAE because this is numerically the  less robust for initialization and re-initialization.\nnc::Int: Number of constraint equations (= length(fc)).\n If structureOfDAE = ImplicitIndexOneDAE, then  nc > 0 signals that constraint equations are present (but they need not to be   at the end of the residue vector).\n If structureOfDAE = LinearDerivativesWithConstraints, then  nc defines the number of constraint equations and these equations must be  at the end of the residue vector.\n If structureOfDAE = ExplicitDerivativesWithoutConstraints, nc = 0 required.\nnz::Int: Number of event indicators\nnw::Int: Number of auxiliary variables (Float64 variables that are additionally computed            and stored at communication points, and where start values can be provided            for initialization)\nzDir::Vector{Int}: Interpretation of event indictors:  zDir[i] = 0: Root is reported for both crossing directions,          = 1: Root is reported when crossing from negative to positive direction          = -1: Root is reported when crossing from positive to negative direction\nx_fixed::Vector{Bool}: = true, x_start[i] is fixed during initialization. = false, x_start[i] might be changed, e.g., due to an initial impulse.\nx_nominal::Vector{Float64}: Nominal values of x. if x_nominal[i]=NaN no nominal value is defined for x[i] and a nominal value is computed (x_nominal[i] = abs(x_start[i]) > 1e-7 ? abs(x_start[i]) : 1.0).\nx_errorControl::Vector{Bool}: = true, the absolute error tolerance is set to 0.1 * relativeTolerance * x_nominal[i]. = false, the absolute error tolerance is switched off (is set to a large value). This is recommended for variables that are basically not limited (for example the angle of a shaft that is permantently rotating in the same direction and therefore becomes larger and larger).\nhev::Float64: Stepsize used during initialization and at event restart if structureOfDAE = ModiaMath.ImplicitIndexOneDAE. Otherwise hev is ignored.\nscaleConstraintsAtEvents::Bool: only kept for backwards compatibility (option is ignored).\njac: Sparse Jacobian datastructure (currently not supported).\nmaxSparsity::Float64: A sparse Jacobian is only used during simulation if sparseness of jac < maxSparsity (currently not supported)\ngetResultNames::Function: Function that returns the names of the variables to be stored in the result data structure.\nstoreResult!::Function: Function that stores the raw results.\ngetResult::Function: Function that resturns the result data structure after the simulation.\ndefaultTolerance::Float64: Model specific default relative tolerance, if not redefined in the call to ModiaMath.simulate!.\ndefaultStartTime::Float64: Model specific default start time in [s], if not redefined in the call to ModiaMath.simulate!.\ndefaultStopTime::Float64: Model specific default stop time in [s], if not redefined in the call to ModiaMath.simulate!.\ndefaultInterval::Float64: Model specific default interval in [s], if not redefined in the call to ModiaMath.simulate!. Result data is stored every defaultInterval seconds. If defaultInterval=NaN, the default interval is computed as interval = (stopTime - startTime)/500.0. \n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.StructureOfDAE",
    "page": "DAE",
    "title": "ModiaMath.DAE.StructureOfDAE",
    "category": "type",
    "text": "@enum StructureOfDAE\n      ImplicitIndexOneDAE\n      LinearDerivativesWithConstraints\n      ExplicitDerivativesWithoutConstraints\n\nEnumeration defining the structure of the DAE of the simulation model.  The following DAE structures are supported  (function getModelResidues!(model, t, x, derx, r, w) returns the residues r):\n\nModiaMath.ImplicitIndexOneDAE\n\nbeginalign\n     z = f_z(xt) \n 0 = r = left beginarrayl\n                    f_d(dotxxtz_i0) \n                    f_c(xtz_i0)\n                  endarray right \n     J = left fracpartial f_dpartial dotx  \n                  fracpartial f_cpartial x right  textis regular (matrix is invertible)\nendalign\n\nInitialization and re-initialization is performed by using an implicit Euler step. When appropriately scaling r, using a step size that tends to zero and under further assumptions, then this solution can be interpreted as analytically integrating over the time instant. This might mean to integrate over Dirac impulses (in case x is discontinuous at this time instant). Since the selected step size is not close to zero, the implicit Euler step will give a very rough approximation of the analytical integral. A much better approximation is achieved with option  LinearDerivativesWithConstraints below where a step size of zero is used.\n\nModiaMath.LinearDerivativesWithConstraints\n\nbeginalign\n     z = f_z(xt) \n 0 = r = left beginarrayl\n                    M_d(xtz_i0) cdot dotx + b_d(xtz_i0)\n                    f_c(xtz_i0)\n                  endarray right \n     J = left M_d  \n                  fracpartial f_cpartial x right  textis regular (matrix is invertible)\nendalign\n\nWhen instantiating a SimulationState, the dimension nc = length(fc) must be provided (nc = 0 is included as a special case).\n\nInitialization and re-initialization is performed by analytically integrating over the initial time or the event time. This might mean to integrate over Dirac impulses (in case x is discontinuous at this time instant). Under certain assumptions a numerical approximation of the mathematically unique solution is computed.\n\nModiaMath.ExplicitDerivativesWithoutConstraints\n\nbeginalign\n     z = f_z(x t) \n 0 = r = f(x t z_i  0) - dotx\nendalign\n\nInitialization and re-initialization is trivial:\n\n     x = x_start \nder(x) = r         # getModelResidues!(model, t, x, zeros(nx), r, w)\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.change!-Tuple{ModiaMath.DAE.SimulationState,Int64,Float64,String}",
    "page": "DAE",
    "title": "ModiaMath.DAE.change!",
    "category": "method",
    "text": "ModiaMath.change!(sim, nr, crossing, crossingAsString; restart = ModiaMath.Restart)\n\nTrigger an event, whenever crossing > 0 changes from false to true or from true to false. The function returns always false.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.defaultVariableName-Tuple{Any,ModiaMath.DAE.VariableCategory,Int64}",
    "page": "DAE",
    "title": "ModiaMath.DAE.defaultVariableName",
    "category": "method",
    "text": "defaultVariableName(model, vcat, vindex)\n\nReturn default names for the variables (e.g. x[1])\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.edge!-Tuple{ModiaMath.DAE.SimulationState,Int64,Float64,String}",
    "page": "DAE",
    "title": "ModiaMath.DAE.edge!",
    "category": "method",
    "text": "ModiaMath.edge!(sim, nr, crossing, crossingAsString; restart = ModiaMath.Restart)\n\nTrigger an event, whenever crossing > 0 switches from false to true. The function returns always false.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.getResultNames-Tuple{Any}",
    "page": "DAE",
    "title": "ModiaMath.DAE.getResultNames",
    "category": "method",
    "text": "getResultNames(model)\n\nReturn a vector of result names.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.getStartTime-Tuple{ModiaMath.DAE.SimulationState}",
    "page": "DAE",
    "title": "ModiaMath.DAE.getStartTime",
    "category": "method",
    "text": "ModiaMath.getStartTime(m::ModiaMath.[AbstractSimulationModel|SimulationState])\n\nReturn startTime of the actual simulation run.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.getStopTime-Tuple{ModiaMath.DAE.SimulationState}",
    "page": "DAE",
    "title": "ModiaMath.DAE.getStopTime",
    "category": "method",
    "text": "ModiaMath.getStopTime(m::ModiaMath.[AbstractSimulationModel|SimulationState])\n\nReturn stopTime of the actual simulation run (when the simulation will be terminated).\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.getTime-Tuple{ModiaMath.DAE.SimulationState}",
    "page": "DAE",
    "title": "ModiaMath.DAE.getTime",
    "category": "method",
    "text": "ModiaMath.getTime(m::ModiaMath.[AbstractSimulationModel|SimulationState])\n\nReturn actual simulation time.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.getTolerance-Tuple{ModiaMath.DAE.SimulationState}",
    "page": "DAE",
    "title": "ModiaMath.DAE.getTolerance",
    "category": "method",
    "text": "ModiaMath.getTolerance(m::ModiaMath.[AbstractSimulationModel|SimulationState])\n\nReturn (relative) tolerance of the actual simulation run.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.getVariableName",
    "page": "DAE",
    "title": "ModiaMath.DAE.getVariableName",
    "category": "function",
    "text": "getVariableName(model, vcat, vindex, nx=0; \n                xNames   =nameVector(\"x\", nx),\n                derxNames=fcNameVector(\"der\", xNames),\n                wNames   =String[])\n\nGiven category vcat and index vindexof category,  return the full path name of the respective variable.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.isAfterSimulationStart-Tuple{ModiaMath.DAE.SimulationState}",
    "page": "DAE",
    "title": "ModiaMath.DAE.isAfterSimulationStart",
    "category": "method",
    "text": "ModiaMath.isAfterSimulationStart(m::ModiaMath.[AbstractSimulationModel|SimulationState])\n\nReturn true, if after start of simulation and false if during initialization.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.isEvent-Tuple{ModiaMath.DAE.SimulationState}",
    "page": "DAE",
    "title": "ModiaMath.DAE.isEvent",
    "category": "method",
    "text": "ModiaMath.isEvent(m::ModiaMath.[AbstractSimulationModel|SimulationState])\n\nReturn true, if event phase of simulation (including initialization)\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.isInitial-Tuple{ModiaMath.DAE.SimulationState}",
    "page": "DAE",
    "title": "ModiaMath.DAE.isInitial",
    "category": "method",
    "text": "ModiaMath.isInitial(m::ModiaMath.[AbstractSimulationModel|SimulationState])\n\nReturn true, if initialization phase of simulation.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.isStoreResult-Tuple{ModiaMath.DAE.SimulationState}",
    "page": "DAE",
    "title": "ModiaMath.DAE.isStoreResult",
    "category": "method",
    "text": "ModiaMath.isStoreResult(m::ModiaMath.[AbstractSimulationModel|SimulationState])\n\nReturn true, if communication point and variables shall be stored in the result data structure.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.isTerminal-Tuple{ModiaMath.DAE.SimulationState}",
    "page": "DAE",
    "title": "ModiaMath.DAE.isTerminal",
    "category": "method",
    "text": "ModiaMath.isTerminal(m::ModiaMath.[AbstractSimulationModel|SimulationState])\n\nReturn true, if terminal phase of simulation.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.isZeroCrossing-Tuple{ModiaMath.DAE.SimulationState}",
    "page": "DAE",
    "title": "ModiaMath.DAE.isZeroCrossing",
    "category": "method",
    "text": "ModiaMath.isZeroCrossing(m::ModiaMath.[AbstractSimulationModel|SimulationState])\n\nReturn true, if simulation model shall compute zero-crossing functions (as required by the integrator).\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.negative!-Tuple{ModiaMath.DAE.SimulationState,Int64,Float64,String}",
    "page": "DAE",
    "title": "ModiaMath.DAE.negative!",
    "category": "method",
    "text": "ModiaMath.negative!(sim, nr, crossing, crossingAsString; restart = ModiaMath.Restart)\n\nReturn crossing < 0 such that a state event is triggered whenever crossing < 0 changes its value.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.positive!-Tuple{ModiaMath.DAE.SimulationState,Int64,Float64,String}",
    "page": "DAE",
    "title": "ModiaMath.DAE.positive!",
    "category": "method",
    "text": "ModiaMath.positive!(sim, nr, crossing, crossingAsString; restart = ModiaMath.Restart)\n\nReturn crossing > 0 such that a state event is triggered whenever crossing > 0 changes its value.\n\nNote, a small hysteresis is added to the crossing function during continuous integration,  in order to reduce issues during event iterations due to small numerical errors. However, at an event instant, crossing > 0 is returned without any hysteresis.\n\nArguments\n\nsim::ModiaMath.SimulationState: (Internal) simulation state provided by simulation engine.\nnr::Int: (> 0 required) Every call of positive!(..) must be identified by a unique value of nr       This value is used as index in a vector that holds the internal memory for crossing > 0.\ncrossing::Float64: Zero crossing function.\ncrossingAsString::String: crossing as string representation. This string is used for log messages.\nrestart::ModiaMath.EventRestart: Restart action after the crossing > 0 event occurred.\n\nExample\n\nfunction computeVariables!(m::Model, sim::SimulationState)\n   ...\n   # f = s > 0 ? fMax : 0.0\n   m.sPos.value = ModiaMath.positive!(sim, 1, m.s.value, \"s\")\n   m.f.value    = m.sPos.value ? m.fMax.value : 0.0\n   ...\nend\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.setNextEvent!-Tuple{ModiaMath.DAE.SimulationState,Float64}",
    "page": "DAE",
    "title": "ModiaMath.DAE.setNextEvent!",
    "category": "method",
    "text": "ModiaMath.setNextEvent!(sim, nextEventTime; integratoToEvent=true, restart = ModiaMath.Restart)\n\nAt an event instant (isEvent(sim) == true) trigger the next time event. \n\nArguments\n\nsim::ModiaMath.SimulationState: (Internal) simulation state provided by simulation engine.\nnextEventTime::Float64: Time instant of the next time event.\nintegrateToEvent::Bool: If true, the integrator integrates exactly to this event. If false, the integrator might integrate beyond nextEventTime (so the step size of the integrator is not influenced by this event). This option is only useful, if information is inquired at the event and restart=ModiaMath.NoRestart is used. The benefit is that the integrator is  practically not influenced by this event.\nrestart::ModiaMath.EventRestart: Restart action after the event at which setNextEvent! was called.\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#ModiaMath.DAE.setNominal!-Tuple{ModiaMath.DAE.SimulationState,Array{Float64,1}}",
    "page": "DAE",
    "title": "ModiaMath.DAE.setNominal!",
    "category": "method",
    "text": "ModiaMath.setNominal!(sim::ModiaMath.SimulationState, x_nominal::Vector{Float64})\n\nAt initialization provide x_nominal, the nominal values of the x-vector, and  store them insim`, the simulation state. These values are used to compute the absolute tolerances of the x-vector for the integrator.\n\nExample\n\nfunction getModelResidues!(m::Model, t, x, derx, r, w)\n   sim = m.simulationState\n   if ModiaMath.isInitial(sim)\n      ModiaMath.setNominal!(sim, [1.0, 1e-5, 1e8])   # if length(x)=3\n   end\n   ...\nend\n\n\n\n\n\n"
},

{
    "location": "lib/DAE.html#DAE-1",
    "page": "DAE",
    "title": "DAE",
    "category": "section",
    "text": "Modules = [ModiaMath.DAE]\r\nPrivate = false"
},

{
    "location": "lib/NonlinearEquations.html#",
    "page": "NonlinearEquations",
    "title": "NonlinearEquations",
    "category": "page",
    "text": ""
},

{
    "location": "lib/NonlinearEquations.html#ModiaMath.NonlinearEquations",
    "page": "NonlinearEquations",
    "title": "ModiaMath.NonlinearEquations",
    "category": "module",
    "text": "module ModiaMath.NonlinearEquations\n\nThis module contains functions to solve nonlinear algebraic equations:\n\nsolveOneNonlinearEquation: Function that computes the solution of one non-linear algebraic equation y=f(u) in one unknown u in a reliable and efficient way using Brents algorithm.\nKINSOL: Module containing functions to solve a system of nonlinear systems of equations with Sundials KINSOL. The module is designed so that the same system is solved several times with KINSOL as needed by Modia simulations (auxiliary memory is only allocated once and not for every call). KINSOL is used in ModiaMath to solve nonlinear algebraic equations during initialization and at events of a simulation.\n\nMain developer\n\nMartin Otter, DLR - Institute of System Dynamics and Control\n\n\n\n\n\n"
},

{
    "location": "lib/NonlinearEquations.html#ModiaMath.NonlinearEquations.solveOneNonlinearEquation-Tuple{Function,Any,Any}",
    "page": "NonlinearEquations",
    "title": "ModiaMath.NonlinearEquations.solveOneNonlinearEquation",
    "category": "method",
    "text": "y = solveOneNonlinearEquation(f, u_min, u_max; u_nominal=1.0, tolerance=100.0*eps())\n\nThis function determines the solution of one non-linear algebraic equation y=f(u) in one unknown u in a reliable and efficient way. It is one of the best numerical algorithms for this purpose. As input, the nonlinear function f(u) has to be given, as well as an interval u_min, u_max that contains the solution, that is f(u_min) and f(u_max) must have a different sign.  The function computes iteratively smaller intervals in which a sign change is present and terminates when the following test is fulfilled (u_nominal is the nominal value of u, so the order of magnitude of u around the solution and tolerance is the relative tolerance; for example, tolerance=1e-10 means that the solution is correct up to 10 digits):\n\nabsTol = 0.1*u_nominal*tolerance\nif abs(length_of_u_interval) <= max(tolerance*abs(u), absTol) || f(u) == 0.0\n    # root u found (interval is small enough)\n\nThe interval reduction is performed using inverse quadratic interpolation (interpolating with a quadratic polynomial through the last three points and computing the zero). If this fails, bisection is used, which always reduces the interval by a factor of two. The inverse quadratic interpolation method has superlinear convergence. This is roughly the same convergence rate as a globally convergent Newton method, but without the need to compute derivatives of the non-linear function. The solver function is a direct mapping of the Algol 60 procedure zero to Julia, from:\n\nBrent R.P. (1973): Algorithms for Minimization without derivatives. Prentice Hall, 1973, pp. 58-59.\nDownload: http://wwwmaths.anu.edu.au/~brent/pd/rpb011i.pdf\nErrata and new print: http://wwwmaths.anu.edu.au/~brent/pub/pub011.html\n\nExamples\n\nimport ModiaMath: solveOneNonlinearEquation\n\nfun(u; w=1.0) = 3*u - sin(w*u) - 1\nu_zero = solveOneNonlinearEquation(u->fun(u; w=3.0), 0.0, 5.0)\nprintln(\"residue = \", fun(u_zero; w=3.0))\n\nRemarks\n\nThe interface was made such that it is identical to function Modelica.Math.Nonlinear.solveOneNonlinearEquation in order that automatic translation of Modelica to Modia is simplified. However, the termination condition was changed: The original Brent algorithm uses an absolute tolerance for termination (abs(lengthofinterval) <= 2eps()abs(u) + tolerance || f(u) == 0.0), whereas this was changed here to take a relative tolerance into account and use a similar definition of tolerances as used by Modia and Modelica, because easier to understand by a user (with relative tolerance and nominal value, instead of relative and absolute tolerances).\nNewer algorithms for the problem are presented in Alefeld, Potra, Shi (1995): Algorithm 748: enclosing zeros of continuous functions. Here, an inverse cubic interpolation is used instead of an inverse quadratic interpolation as in case of the Brent algorithm. The numerical experiments in this article with 15 test problems show that Brents algorithm needs about 5% more function evaluations for all the 15 test problems as the newly presented algorithm 4.2 (in some cases Brents algorithm needs slightly less and in other cases slightly more function evaluations).\nJulia package Roots.jl provides various algorithms to solve the problem above but it seems that Brents algorithm is not yet included (as of Nov. 2018).\n\n```\n\n\n\n\n\n"
},

{
    "location": "lib/NonlinearEquations.html#ModiaMath.NonlinearEquations.KINSOL",
    "page": "NonlinearEquations",
    "title": "ModiaMath.NonlinearEquations.KINSOL",
    "category": "module",
    "text": "module KINSOL - Solve nonlinear equation system with Sundials KINSOL\n\nThe goal is to solve the same system several times with KINSOL.\n\n\n\n\n\n"
},

{
    "location": "lib/NonlinearEquations.html#NonlinearEquations-1",
    "page": "NonlinearEquations",
    "title": "NonlinearEquations",
    "category": "section",
    "text": "Modules = [ModiaMath.NonlinearEquations, ModiaMath.NonlinearEquations.KINSOL]\r\nPrivate = false"
},

{
    "location": "lib/Variables.html#",
    "page": "Variables",
    "title": "Variables",
    "category": "page",
    "text": ""
},

{
    "location": "lib/Variables.html#ModiaMath.Variables",
    "page": "Variables",
    "title": "ModiaMath.Variables",
    "category": "module",
    "text": "module ModiaMath.Variables\n\nProvide variable types for simulation models.\n\nThis module provides functions to declare simulation Variables as special structs that have a value and associated attributes. Furthermore, there are functions to copy the integrator interface variables (x, derx, residue) to the appropriate Variables and  vice versa. Therefore, the modeler does not have to care about, how the values of the variables are mapped to the integrator interface. Typically, a model is constructed with macro @component using RealXXX variable declarations. Example:\n\n  using ModiaMath\n  using StaticArrays\n\n  @component Pendulum(;L=1.0, m=1.0, d=0.1, g=9.81) begin\n     phi = RealScalar(start=pi/2, unit=\"rad\"    , fixed=true,               numericType=ModiaMath.XD_EXP)\n     w   = RealScalar(start=0.0 , unit=\"rad/s\"  , fixed=true, integral=phi, numericType=ModiaMath.XD_EXP)\n     a   = RealScalar(            unit=\"rad/s^2\",             integral=w  , numericType=ModiaMath.DER_XD_EXP) \n     r   = RealSVector{2}(        unit=\"m\"      ,                           numericType=ModiaMath.WC)\n  end;\n\nThe following variable types are currently supported:\n\nVariable types Description\nRealScalar Scalar variable with Float64 value\nRealSVector{Size} Variable with SVector{Size,Float64} value\nRealSVector3 Variable with SVector{3,Float64} value\nRealVariable{VType, EType} Variable with value/element type VType/EType\n\nAll of them have the following attributes:\n\nVariable attributes Description\nvalue::VType Value of the variable (scale, vector, matrix, ...)\ninfo::AbstractString Description text\ncausality::ModiaMath.Causality Causality of variable\nvariability::ModiaMath.Variability Variability of variable\nstart::EType Start value of variable\nfixed::Bool fixed = true: value=start after initialization\nanalysis::VariableAnalysisType Analysis for which the variable is used\nmin::EType Minimum value of value and of start\nmax::EType Maximum value of value and of start\nnominal::EType Nominal value of value (used to improve numerics)\nflow::Bool = true if variable is a flow variable\nnumericyType::ModiaMath.NumericType Defines how variable is used by the integrator\nintegral = the integral variable or nothing\nunit::String Unit of value and of start\n\nThe following functions are provided to perform the actual copy-operations and/or to inquire information about the variables in the model\n\nModiaMath.copy_start_to_x!\nModiaMath.copy_x_and_derx_to_variables!\nModiaMath.copy_variables_to_residue!\nModiaMath.print_ModelVariables\nModiaMath.get_xTable\nModiaMath.get_copyToVariableTable\nModiaMath.get_copyToResidueTable\nModiaMath.get_copyToResultTable\n\nMain developer\n\nMartin Otter,  DLR - Institute of System Dynamics and Control\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.AnalysisType",
    "page": "Variables",
    "title": "ModiaMath.Variables.AnalysisType",
    "category": "type",
    "text": "@enum AnalysisType KinematicAnalysis QuasiStaticAnalysis DynamicAnalysis\n\nType of analyis that is actually carried out. The AnalysisType is set by the user of the simulation model. Variables are declared in the model with VariableAnalysisType. Variables with [VariableAnalysisType](@ref) <= AnalysisTypeare removed from the analysis and do not show up in the result. For example, an *acceleration* would be declared asOnlyDynamicAnalysisand then this variable would not show up in the result, ifAnalysisType = KinematicAnalysisorQuasiStaticAnalysis`.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.Causality",
    "page": "Variables",
    "title": "ModiaMath.Variables.Causality",
    "category": "type",
    "text": "@enum Causality Parameter CalculatedParameter Input Output Local Independent\n\nCausality of variable (only used when connecting models)\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.ModelVariables",
    "page": "Variables",
    "title": "ModiaMath.Variables.ModelVariables",
    "category": "type",
    "text": "vars = ModiaMath.ModelVariables(model::ModiaMath.AbstractComponentWithVariables; \n                                analysis::AnalysisType = ModiaMath.DynamicAnalysis)\n\nReturn a struct that contains all the variables in model in a form so that fast copying from the integrator interface to the variables can be performed, and vice versa.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.NumericType",
    "page": "Variables",
    "title": "ModiaMath.Variables.NumericType",
    "category": "type",
    "text": "@enum NumericType XD_EXP XD_IMP ..\n\nDefines how a variable is used in the integrator interface. The goal is to describe an implicit index-1 DAE:\n\nbeginalign\n z = f_z(x t) \n 0 = f_d(dotx x t z_i  0) \n 0 = f_c(x t z_i  0) \n J = left fracpartial f_dpartial dotx  \n              fracpartial f_cpartial x right  textis regular\nendalign\n\nThe integrator interface to the model is\n\ngetModelResidues(m::AbstractSimulationModel, t::Float64, x::Vector{Float64},  \n                 derx::Vector{Float64}, r::Vector{Float64}\n\nIn the following table it is shown how NumericType values of variables are mapped to this integrator interface:\n\n@enum value Mapping/Description\nXD_EXP Copied from x into variable; der(v) is computed by model\nXD_IMP Copied from x into variable; der(v) copied from derx\nXA Copied from x into variable; derx is ignored (algebraic variable)\nDER_XD_EXP Computed by model; automatic: residue = der(v) - derx\nDER_XD_IMP Copied from derx into variable (= der(v); v has type XD_IMP\nLAMBDA Copied from derx into variable; index-1 algebraic variable\nMUE Copied from derx; stabilizing variable to reduced index to index 1\nFD_IMP Is copied from variable to residue; part of f_d equations\nFC Is copied from variable to residue; part of f_c equations\nWR Computed by model; stored in result; used to compute residues\nWC Computed by model; stored in result; only available at communication points\nTIME Copied from t into variable; independent variable time\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.RealSVector",
    "page": "Variables",
    "title": "ModiaMath.Variables.RealSVector",
    "category": "type",
    "text": "v = RealSVector{Size}(..)\n\nGenerate a variable of type RealVariable{SVector{Size,Float64},Float64} where the values have type StaticArrays.SVector{Size,Float64}. The argument list is described in module Variables.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.RealSVector3",
    "page": "Variables",
    "title": "ModiaMath.Variables.RealSVector3",
    "category": "type",
    "text": "v = RealSVector3(..)\n\nGenerate a variable of type RealVariable{SVector{3,Float64},Float64} where the values have type StaticArrays.SVector{3,Float64}. The argument list is described in module Variables.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.RealScalar",
    "page": "Variables",
    "title": "ModiaMath.Variables.RealScalar",
    "category": "type",
    "text": "v = RealScalar(..)\n\nGenerate a variable of type RealVariable{Float64,Float64} to define scalar, Float64, real variables. The argument list is described in module Variables.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.RealVariable",
    "page": "Variables",
    "title": "ModiaMath.Variables.RealVariable",
    "category": "type",
    "text": "v = RealVariable{ValueType,ElementType}(...)\n\nGenerate a new variable with v.value::ValueType, v.nominal::ElementType. The argument list is described in module Variables.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.SimulationModel",
    "page": "Variables",
    "title": "ModiaMath.Variables.SimulationModel",
    "category": "type",
    "text": "sm = ModiaMath.SimulationModel(model::ModiaMath.AbstractComponentWithVariables;\n                               startTime = 0.0, stopTime = 1.0, tolerance = 1e-4,\n                               interval = (stopTime-startTime)/500.0)\n\nReturn a simulationModel sm struct that can be simulated with ModiaMath.simulate!. The given model is typically constructed with the @component macro. As keyword arguments, default values for startTime, stopTime, tolerance, interval can be given.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.Variability",
    "page": "Variables",
    "title": "ModiaMath.Variables.Variability",
    "category": "type",
    "text": "@enum Variability Constant Fixed Tunable Discrete Continuous\n\nCurrently not used. Will be used in the future to store only the minimum information in the result (store value only, when it can potentially change).\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.VariableAnalysisType",
    "page": "Variables",
    "title": "ModiaMath.Variables.VariableAnalysisType",
    "category": "type",
    "text": "@enum VariableAnalysisType AllAnalysis QuasiStaticAndDynamicAnalysis \n                           OnlyDynamicAnalysis NotUsedInAnalysis\n\nType of analysis that can be carried out with the variable (e.g. a force would be defined as QuasiStaticAndDynamicAnalysis, an acceleration would be defined as OnlyDynamicAnalysis and a position variable would be defined as AllAnalysis.\n\nVariables with VariableAnalysisType <=AnalysisType are removed from the analysis and do not show up in the result. For example, an acceleration would be declared as OnlyDynamicAnalysis and then this variable would not show up in the result, if AnalysisType = KinematicAnalysis or QuasiStaticAnalysis.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.componentName-Tuple{ModiaMath.AbstractComponentWithVariables}",
    "page": "Variables",
    "title": "ModiaMath.Variables.componentName",
    "category": "method",
    "text": "name = ModiaMath.componentName(\n           component::ModiaMath.AbstractComponentWithVariables)\n\nReturn name of component (without the leading path) as Symbol. \n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.copy_start_to_x!-Tuple{ModiaMath.Variables.ModelVariables,Array{Float64,1},Array{Bool,1},Array{Float64,1}}",
    "page": "Variables",
    "title": "ModiaMath.Variables.copy_start_to_x!",
    "category": "method",
    "text": "copy_start_to_x!(vars::ModelVariables, x::Vector{Float64}, x_fixed::Vector{Bool}, x_nominal::Vector{Float64})\n\nCopy start, fixed and nominalvalues of variables vars to x, x_fixed, and x_nominal vectors.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.copy_variables_to_residue!-Tuple{ModiaMath.Variables.ModelVariables,Array{Float64,1},Array{Float64,1},Array{Float64,1}}",
    "page": "Variables",
    "title": "ModiaMath.Variables.copy_variables_to_residue!",
    "category": "method",
    "text": "copy_variables_to_residue!(vars::ModelVariables, x::Vector{Float64},\n                           derx::Vector{Float64}, residue::Vector{Float64})\n\nCopy variables vars to residue vector and include the inherent residue equations of XD_EXP variables (residue = der(v) - derx).\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.copy_x_and_derx_to_variables!-Tuple{Float64,Array{Float64,1},Array{Float64,1},ModiaMath.Variables.ModelVariables}",
    "page": "Variables",
    "title": "ModiaMath.Variables.copy_x_and_derx_to_variables!",
    "category": "method",
    "text": "copy_x_and_derx_to_variables!(time::Float64, x::Vector{Float64}, \n                              derx::Vector{Float64}, vars::ModelVariables)\n\nCopy xand derxof the integrator interface to the model variables vars.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.fullName-Tuple{ModiaMath.AbstractComponentWithVariables}",
    "page": "Variables",
    "title": "ModiaMath.Variables.fullName",
    "category": "method",
    "text": "name = ModiaMath.fullName(component::ModiaMath.AbstractComponentWithVariables)\n\nReturn full path name of component (including root name) as Symbol.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.get_copyToResidueTable-Tuple{ModiaMath.Variables.ModelVariables}",
    "page": "Variables",
    "title": "ModiaMath.Variables.get_copyToResidueTable",
    "category": "method",
    "text": "table = ModiaMath.get_copyToResidueTable(vars::ModelVariables)\n\nFunction returns a DataFrames tables of all the variables in vars that are copied from variables to the residue vector.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.get_copyToResultTable-Tuple{ModiaMath.Variables.ModelVariables}",
    "page": "Variables",
    "title": "ModiaMath.Variables.get_copyToResultTable",
    "category": "method",
    "text": "table = ModiaMath.get_copyToResultTable(vars::ModelVariables)\n\nFunction returns a DataFrames tables of all the variables in vars that are copied from variables to the result.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.get_copyToVariableTable-Tuple{ModiaMath.Variables.ModelVariables}",
    "page": "Variables",
    "title": "ModiaMath.Variables.get_copyToVariableTable",
    "category": "method",
    "text": "table = ModiaMath.get_copyToVariableTable(vars::ModelVariables)\n\nFunction returns a DataFrames tables of all the variables in vars that are copied from x/derx to the variables.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.get_xTable-Tuple{ModiaMath.Variables.ModelVariables}",
    "page": "Variables",
    "title": "ModiaMath.Variables.get_xTable",
    "category": "method",
    "text": "table = ModiaMath.get_xTable(vars::ModelVariables)\n\nFunction returns a DataFrames tables of all the variables stored in x-vector in vars.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.initComponent!-Tuple{ModiaMath.AbstractComponentWithVariables,Any,Any}",
    "page": "Variables",
    "title": "ModiaMath.Variables.initComponent!",
    "category": "method",
    "text": "initComponent!(within::ModiaMath.AbstractComponentWithVariables, component, name)\n\nInitialize component with name (of String or Symbol type) for component within. If component::ModiaMath.AbstractComponentWithVariables the within object and name are stored in component._internal.\n\nAdditionally, and for all other types of component the following statement is executed:\n\n   setfield!(within, Symbol(name), component)\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.instanceName-Tuple{ModiaMath.AbstractComponentWithVariables}",
    "page": "Variables",
    "title": "ModiaMath.Variables.instanceName",
    "category": "method",
    "text": "ModiaMath.instanceName(component::ModiaMath.AbstractComponentWithVariables)\n\nReturn instance name of component (= full name but without root name) as Symbol.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.print_ModelVariables-Tuple{ModiaMath.Variables.ModelVariables}",
    "page": "Variables",
    "title": "ModiaMath.Variables.print_ModelVariables",
    "category": "method",
    "text": "table = print_ModelVariables(obj)\n\nPrint all the variables in obj in form of DataFrames tables. obj can be of type ModiaMath.ModelVariables or ModiaMath.SimulationModel.\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#ModiaMath.Variables.@component-Tuple{Any,Any}",
    "page": "Variables",
    "title": "ModiaMath.Variables.@component",
    "category": "macro",
    "text": "@component ComponentName(arguments) begin ... end\n\nMacro to generate a mutable struct ComponentName. An instance of this struct can be used as simulationModel in constructor ModiaMath.SimulationModel and can then be simulated with ModiaMath.simulate!.\n\nThe arguments must be keyword arguments and are used as keyword arguments for the generated constructor function ComponentName. The code inbegin ... endis  basically the body of the generated constructor function. All left-hand-side (scalar or vector) symbols present betweenbegin ... end, as well as all keyword-argumentsargumentsare declared as fields in structComponentName`.\n\nExamples\n\nusing ModiaMath\n\n@component Pendulum(;L=1.0, m=1.0, d=0.1, g=9.81, phi0_deg=90.0) begin\n   @assert(L > 0.0)\n   @assert(m > 0.0)\n   @assert(d >= 0.0)\n  \n   phi = RealScalar(..)\n   w   = RealScalar(..)\n   a   = RealScalar(..) \n   r   = RealSVector{2}(..)\nend\n\npendulum = Pendulum(L=1.8)\n\nThe symbols  L, m, d, g, phi0_deg, phi, w, a, r are used as fields in pendulum (so pendulum.L = 1.8).\n\n\n\n\n\n"
},

{
    "location": "lib/Variables.html#Variables-1",
    "page": "Variables",
    "title": "Variables",
    "category": "section",
    "text": "Modules = [ModiaMath.Variables]\r\nPrivate = false"
},

{
    "location": "lib/Logging.html#",
    "page": "Logging",
    "title": "Logging",
    "category": "page",
    "text": ""
},

{
    "location": "lib/Logging.html#ModiaMath.Logging",
    "page": "Logging",
    "title": "ModiaMath.Logging",
    "category": "module",
    "text": "module ModiaMath.Logger\n\nLog model evaluations.\n\nMain developer\n\nMartin Otter, DLR - Institute of System Dynamics and Control\n\n\n\n\n\n"
},

{
    "location": "lib/Logging.html#ModiaMath.Logging.Logger",
    "page": "Logging",
    "title": "ModiaMath.Logging.Logger",
    "category": "type",
    "text": "mutable struct Logger - Log model evaluations\n\n\n\n\n\n"
},

{
    "location": "lib/Logging.html#ModiaMath.Logging.SimulationStatistics",
    "page": "Logging",
    "title": "ModiaMath.Logging.SimulationStatistics",
    "category": "type",
    "text": "mutable struct SimulationStatistics - Collect statistics of the last simulation run.\n\nThe following data is stored in this structure:\n\nstructureOfDAE: Structure of DAE\ncpuTimeInitialization: CPU-time for initialization\ncpuTimeIntegration: CPU-time for integration\nstartTime: start time of the integration\nstopTime: stop time of the integration\ninterval: communication interval of the integration\ntolerance: relative tolerance used for the integration\nnEquations: number of equations (length of y and of yp)\nnConstraints: number of constraint equations\nnResults: number of time points stored in result data structure\nnSteps: number of (successful) steps\nnResidues: number of calls to compute residues (includes residue calls for Jacobian)\nnZeroCrossing: number of calls to compute zero crossings\nnJac: number of calls to compute Jacobian\nnTimeEvents: number of time events\nnRestartEvents: number of events with integrator restart\nnErrTestFails: number of fails of error tests\nh0: stepsize used at the first step\nhMin: minimum integration stepsize\nhMax: maximum integration stepsize\norderMax: maximum integration order\nsparseSolver = true: if sparse solver used, otherwise dense solver\nnGroups: if sparseSolver, number of column groups to compute Jacobian (e.g. if nEquations=100, nGroups=5, then 5+1=6 model evaluations are needed to compute the Jacobian numerically, instead of 101 model evaluations without taking the sparseness structure into account).\n\n\n\n\n\n"
},

{
    "location": "lib/Logging.html#ModiaMath.Logging.isLogEvents-Tuple{ModiaMath.Logging.Logger}",
    "page": "Logging",
    "title": "ModiaMath.Logging.isLogEvents",
    "category": "method",
    "text": "ModiaMath.isLogEvents(obj)\n\nReturn true, if logger settings require to print event messages of the model (obj must be of type ModiaMath.SimulationState, ModiaMath.AbstractSimulationModel,  or ModiaMath.Logger).\n\n\n\n\n\n"
},

{
    "location": "lib/Logging.html#ModiaMath.Logging.isLogInfos-Tuple{ModiaMath.Logging.Logger}",
    "page": "Logging",
    "title": "ModiaMath.Logging.isLogInfos",
    "category": "method",
    "text": "ModiaMath.isLogInfos(obj)\n\nReturn true, if logger settings require to print info messages of the model (obj must be of type ModiaMath.SimulationState, ModiaMath.AbstractSimulationModel,  or ModiaMath.Logger).\n\n\n\n\n\n"
},

{
    "location": "lib/Logging.html#ModiaMath.Logging.isLogProgress-Tuple{ModiaMath.Logging.Logger}",
    "page": "Logging",
    "title": "ModiaMath.Logging.isLogProgress",
    "category": "method",
    "text": "ModiaMath.isLogProgress(logger::ModiaMath.Logger)\n\nReturn true, if logger settings require to print progress messages of the model\n\n\n\n\n\n"
},

{
    "location": "lib/Logging.html#ModiaMath.Logging.isLogStatistics-Tuple{ModiaMath.Logging.Logger}",
    "page": "Logging",
    "title": "ModiaMath.Logging.isLogStatistics",
    "category": "method",
    "text": "ModiaMath.isLogStatistics(logger::ModiaMath.Logger)\n\nReturn true, if logger settings require to print statistics messages of the model.\n\n\n\n\n\n"
},

{
    "location": "lib/Logging.html#ModiaMath.Logging.isLogWarnings-Tuple{ModiaMath.Logging.Logger}",
    "page": "Logging",
    "title": "ModiaMath.Logging.isLogWarnings",
    "category": "method",
    "text": "ModiaMath.isLogWarnings(obj)\n\nReturn true, if logger settings require to print warning messages of the model (obj must be of type ModiaMath.SimulationState, ModiaMath.AbstractSimulationModel,  or ModiaMath.Logger).\n\n\n\n\n\n"
},

{
    "location": "lib/Logging.html#ModiaMath.Logging.logOff!-Tuple{ModiaMath.Logging.Logger}",
    "page": "Logging",
    "title": "ModiaMath.Logging.logOff!",
    "category": "method",
    "text": "ModiaMath.logOff!(obj)\n\nDisable logging on obj (of type ModiaMath.SimulationState, ModiaMath.AbstractSimulationModel,  or ModiaMath.Logger)\n\n\n\n\n\n"
},

{
    "location": "lib/Logging.html#ModiaMath.Logging.logOn!-Tuple{ModiaMath.Logging.Logger}",
    "page": "Logging",
    "title": "ModiaMath.Logging.logOn!",
    "category": "method",
    "text": "ModiaMath.logOn!(obj)\n\nEnable logging on obj (of type ModiaMath.SimulationState, ModiaMath.AbstractSimulationModel,  or ModiaMath.Logger)\n\n\n\n\n\n"
},

{
    "location": "lib/Logging.html#ModiaMath.Logging.setLogCategories!-Tuple{ModiaMath.Logging.Logger,Array{Symbol,1}}",
    "page": "Logging",
    "title": "ModiaMath.Logging.setLogCategories!",
    "category": "method",
    "text": "ModiaMath.setLogCategories!(obj, categories; reinit=true)\n\nSet log categories on obj (of type ModiaMath.SimulationState, ModiaMath.AbstractSimulationModel,  or ModiaMath.Logger) as vector of symbols, e.g. setLogCategories!(simulationModel, [:LogProgess]). Supported categories:\n\n:LogStatistics, print statistics information at end of simulation\n:LogProgress, print progress information during simulation\n:LogInfos, print information messages of the model\n:LogWarnings, print warning messages of the model\n:LogEvents, log events of the model\n\nIf option reinit=true, all previous set categories are reinitialized to be no longer present. If reinit=false, previously set categories are not changed.\n\n\n\n\n\n"
},

{
    "location": "lib/Logging.html#Logging-1",
    "page": "Logging",
    "title": "Logging",
    "category": "section",
    "text": "Modules = [ModiaMath.Logging]"
},

{
    "location": "lib/Frames.html#",
    "page": "Frames",
    "title": "Frames",
    "category": "page",
    "text": ""
},

{
    "location": "lib/Frames.html#ModiaMath.Frames",
    "page": "Frames",
    "title": "ModiaMath.Frames",
    "category": "module",
    "text": "module ModiaMath.Frames\n\nThis module contains functions for frames that is coordinate systems in 3D. The orientation of a frame is described either with a 3x3 rotation matrix  or with a quaternion vector and its origin is described with a Vector3D:\n\nconstModiaMath.RotationMatrix = SMatrix{3,3,Float64,9}`:  Type of a Rotation matrix to rotate from a frame 1 into a frame 2.\nconstModiaMath.Quaternion = SVector{4,Float64}`: Type of a Quaternion vector to rotate from a frame 1 into a frame 2.\nconstModiaMath.Vector3D = SVector{3,Float64}`: Type of a vector in 3D space (e.g. position vector of the origin of a frame).\n\nThe following constants are defined\n\nconstModiaMath.NullRotation:  RotationMatrix with no rotation from a frame 1 into a frame 2.\nconstModiaMath.NullQuaternion: Quaternion vector with no rotation from a frame 1 into a frame 2.\nconstModiaMath.ZeroVector3D: Vector3D with only zero elements.\n\nIf an angle is given as an argument to one of the functions below, it might be a number (interpreted as having unit rad) or a number with a unit (for example: using Unitful; angle = 90u\"°\").\n\nConstructors for a RotationMatrix R\n\nThe following functions return a ModiaMath.RotationMatrixR to rotate a frame 1 into a frame 2.\n\nFunction Description\nModiaMath.rot1(angle) Rotate around angle along x-axis\nModiaMath.rot2(angle) Rotate around angle along y-axis\nModiaMath.rot3(angle) Rotate around angle along z-axis\nModiaMath.rot123(angle1, angle2, angle3) Rotate around angles along x,y,z-axes\nModiaMath.rot_e(e, angle) Rotate around angle along unit vector e\nModiaMath.rot_nxy(nx, ny) nx/ny are in x/y-direction of frame 2\nModiaMath.from_q(q) Return R from Quaternion q\n\nConstructors for a Quaternion q\n\nThe following functions return a ModiaMath.Quaternionq to rotate a frame 1 into a frame 2. Since q and -q define the same rotation the constructor functions have a keyword argument q_guess::Quaternion = NullQuaternion. From the two possible solutions q the one is returned that is closer to q_guess.\n\nFunction Description\nModiaMath.qrot1(angle) Rotate around angle along x-axis\nModiaMath.qrot2(angle) Rotate around angle along y-axis\nModiaMath.qrot3(angle) Rotate around angle along z-axis\nModiaMath.qrot123(angle1, angle2, angle3) Rotate around angles along x,y,z-axes\nModiaMath.qrot_e(e, angle) Rotate around angle along unit vector e\nModiaMath.qrot_nxy(nx, ny) nx/ny are in x/y-direction of frame 2\nModiaMath.from_R(R) Return q from RotationMatrix R\n\nOperations on Frames\n\nThe following functions provide operations on frames. The orientation of a frame is  defined with argument Rq meaning it can be either a  ModiaMath.RotationMatrix R or a  ModiaMath.Quaternion q (to rotate a frame 1 into a frame 2).\n\nFunction Description\nModiaMath.resolve1(Rq, v2) Transform vector v from frame 2 to frame 1\nModiaMath.resolve2(Rq, v1) Transform vector v from frame 1 to frame 2\nModiaMath.absoluteRotation(Rq01, Rq12) Return rotation 0->2 from rot. 0->1 and 1->2\nModiaMath.relativeRotation(Rq01, Rq02) Return rotation 1->2 from rot. 0->1 and 0->2\nModiaMath.inverseRotation(Rq01) Return rotation 1->0 from rot, 0->1\nModiaMath.planarRotationAngle(e,v1,v2) Return angle of planar rotation along e\nModiaMath.eAxis(axis) Return unit vector e in direction of axis\nModiaMath.skew(v) Return skew-symmetric matrix of vector v\n\nInterpolation of Frames\n\nGiven a set of frames by a vector r of position vectors to their origins and  and an optional vector q of Quaternions of their absolute orientations, then the following functions interpolate linearly in these frames:\n\nFunction Description\nModiaMath.Path(r,q) Return path defined by a vector of frames\nModiaMath.t_pathEnd(path) Return path parameter t_end of last frame\nModiaMath.interpolate(path,t) Return (rt,qt) of Path at path parameter t\nModiaMath.interpolate_r(path,t) Return rt of Path at path parameter t\n\nExamples\n\nusing ModiaMath\nusing Unitful\n\n# R1,R2,R3 are the same RotationMatrices\nR1 = ModiaMath.rot1(pi/2)\nR2 = ModiaMath.rot1(90u\"°\")\nR3 = ModiaMath.rot_e([1,0,0], pi/2)\n\nMain developer\n\nMartin Otter,  DLR - Institute of System Dynamics and Control\n\nThe functions of this module are mostly a mapping of some of the functions of the Modelica Standard Library from Modelica  (Modelica.Mechanics.MultiBody.Frames) to Julia (taking advantage of Julia features  such as multiple dispatch and unit package Unitful).\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.NullQuaternion",
    "page": "Frames",
    "title": "ModiaMath.Frames.NullQuaternion",
    "category": "constant",
    "text": "const ModiaMath.NullQuaternion = Quaternion(0,0,0,1)\n\nConstant Quaternion vector of a null rotation (= no rotation from frame 1 to frame 2)\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.NullRotation",
    "page": "Frames",
    "title": "ModiaMath.Frames.NullRotation",
    "category": "constant",
    "text": "Constant RotationMatrix that defines no rotation from frame 1 to frame 2.\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.ZeroVector3D",
    "page": "Frames",
    "title": "ModiaMath.Frames.ZeroVector3D",
    "category": "constant",
    "text": "const ModiaMath.ZeroVector3D = Vector3D(0.0, 0.0, 0.0)\n\nConstant of a Vector3D where all elements are zero\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.Path",
    "page": "Frames",
    "title": "ModiaMath.Frames.Path",
    "category": "type",
    "text": "path = ModiaMath.Path(r::Vector{Vector3D},\n                      q::Vector{Quaternion} = Quaternion[];\n                      v = ones(length(r)))\n\nReturn an instance of a new Path object. The Path object consists of n frames defined by the position vectors of their origins (r[i] for frame i)  and optionally of their absolute rotation quaternions (q[i] for frame i) describing the rotation from the world frame to the respective frame.\n\nA path parameter t is defined in the following way on these frames:\n\nt[1] = 0.0.\nt[i] = t[i-1] + pathLength_i/(v[i]+v[i-1])/2 if the origins of frames i-1 and i do not coincide.\nt[i] = t[i-1] + pathAngle_i/(v[i]+v[i-1])/2 if the origins of frames i-1 and i do coincide.\n\nHereby pathLength_i is the distance between the origins of frames i-1 and i in [m] and pathAngle_i is the planar rotation angle between frames i-1 and i in [rad].\n\nIf v[i] is the desired velocity or angular velocity at frame i, then path parameter t is approximately the time to move along the path. The time instant t_end of the last frame can be inquired with ModiaMath.t_pathEnd(path). For example, if a simulation shall be performed in such  a way that the simulation should start with the first frame and end at stopTime at the last frame, then the path parameter should be selected as t = time*t_end/stopTime.\n\nGiven the actual path parameter, typically 0 <= t <= t_end (if t is outside of this interval, the frame at t is determined by extrapolation through the first two or the last two frames),  the corresponding frame is determined by linear interpolation in the following way:\n\n(rt, qt) = interpolate(  path,t)\n rt      = interpolate_r(path,t)\n\nwhere rt is the position vector to the origin of the frame at path parameter t and qt is the absolute quaternion of the frame at path parameter t.\n\nExample\n\nimport ModiaMath\nusing Unitful\n\nr = [ ModiaMath.Vector3D(1,0,0),\n      ModiaMath.Vector3D(0,1,0),\n      ModiaMath.Vector3D(0,0,1) ]\nq = [ ModiaMath.NullQuaternion,\n      ModiaMath.qrot1(45u\"°\"),\n      ModiaMath.qrot2(60u\"°\")]\n\npath     = ModiaMath.Path(r,q)\nt_end    = ModiaMath.t_pathEnd(path)\ndt       = 0.1\nstopTime = 2.0 \ntime     = 0.0\n\nwhile time <= stopTime\n   (rt, qt) = ModiaMath.interpolate(path, time*t_end/stopTime)\n   time += dt\nend\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.Quaternion",
    "page": "Frames",
    "title": "ModiaMath.Frames.Quaternion",
    "category": "type",
    "text": "const ModiaMath.Quaternion = SVector{4,Float64}\n\nDescribes the rotation from a frame 1 into a frame 2 with a quaternion vector. If e is the (normalized) axis of rotation to rotate frame 1 into frame 2 (either resolved in frame 1 or frame 2) and angle is the rotation angle for this rotation then the quaternion vector q::ModiaMath.Quaternions is defined as:\n\nq = [e*sin(angle/2),\n       cos(angle/2]\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.RotationMatrix",
    "page": "Frames",
    "title": "ModiaMath.Frames.RotationMatrix",
    "category": "type",
    "text": "const ModiaMath.RotationMatrix = SMatrix{3,3,Float64,9}\n\nDescribes the rotation from a frame 1 into a frame 2. An instance R of RotationMatrix has the following interpretation:\n\nR::RotationMatrix = [ex ey ez]   \n\nwhere ex, ey, ez are unit vectors in the direction of the x-axis, y-axis, and z-axis of frame 1, resolved in frame 2, respectively (for example ex=[1.0, 0.0, 0.0]) Therefore, if v1 is vector v resolved in frame 1 and v2 is vector v resolved in frame 2, the following relationship holds: \n\nv2 = R*v1   \nv1 = R\'*v2\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.Vector3D",
    "page": "Frames",
    "title": "ModiaMath.Frames.Vector3D",
    "category": "type",
    "text": "const ModiaMath.Vector3D = SVector{3,Float64}\n\nType of a vector in 3D space (e.g. position vector of the origin of a frame)\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.absoluteRotation-Tuple{StaticArrays.SArray{Tuple{3,3},Float64,2,9},StaticArrays.SArray{Tuple{3,3},Float64,2,9}}",
    "page": "Frames",
    "title": "ModiaMath.Frames.absoluteRotation",
    "category": "method",
    "text": " R2 = ModiaMath.absoluteRotation(R1, R_rel) \n q2 = ModiaMath.absoluteRotation(q1, q_rel)\n\nReturn ModiaMath.RotationMatrixR2 or ModiaMath.Quaternionq2 defining the rotation from frame 0 to frame 2 from RotationMatrix R1 or Quaternion q1that define the rotation from frame 0 to frame 1 and the relative RotationMatrix R_rel or the relative Quaternion q_rel that define the rotation from frame 1 to frame 2.\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.assertQuaternion-Tuple{AbstractArray{T,1} where T}",
    "page": "Frames",
    "title": "ModiaMath.Frames.assertQuaternion",
    "category": "method",
    "text": "ModiaMath.assertQuaternion(q::AbstractVector)\n\nAssert that vector q has the properties of a Quaternion vector (has 4 elements, norm(q) = 1)\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.assertRotationMatrix-Tuple{AbstractArray{T,2} where T}",
    "page": "Frames",
    "title": "ModiaMath.Frames.assertRotationMatrix",
    "category": "method",
    "text": "ModiaMath.assertRotationMatrix(R::AbstractMatrix)\n\nAssert that matrix R has the properties of a rotation matrix (is 3x3 and R\'*R - eye(3) = zeros(3,3))\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.eAxis-Tuple{Int64}",
    "page": "Frames",
    "title": "ModiaMath.Frames.eAxis",
    "category": "method",
    "text": "e = eAxis(axis::Int)\n\nReturn unit vector e in direction of axis axis (axis = 1,2,3 or -1,-2-,3).\n\nExample\n\nimport ModiaMath\n\ne1 = ModiMath.eAxis(1)    # e1 = Vector3D(1.0,  0.0, 0.0)\ne2 = ModiMath.eAxis(-2)   # d2 = Vector3D(0.0, -1.0, 0.0)\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.from_R-Tuple{StaticArrays.SArray{Tuple{3,3},Float64,2,9}}",
    "page": "Frames",
    "title": "ModiaMath.Frames.from_R",
    "category": "method",
    "text": "q = ModiaMath.from_R(R::ModiaMath.RotationMatrix;\n                     q_guess = NullQuaternion)\n\nReturn Quaternion q from RotationMatrix R.\n\nFrom the two possible solutions q the one is returned that is closer  to q_guess (note, q and -q define the same rotation).\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.from_q-Tuple{StaticArrays.SArray{Tuple{4},Float64,1,4}}",
    "page": "Frames",
    "title": "ModiaMath.Frames.from_q",
    "category": "method",
    "text": "R = ModiaMath.from_q(q::ModiaMath.Quaternion)\n\nReturn RotationMatrix R from Quaternion q.\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.interpolate-Tuple{ModiaMath.Frames.Path,Number}",
    "page": "Frames",
    "title": "ModiaMath.Frames.interpolate",
    "category": "method",
    "text": "(rt, qt) = ModiaMath.interpolate(path, t)\n\nReturn position rtand Quaternion qt of path::ModiaMath.Path at path parameter t::Number.\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.interpolate_r-Tuple{ModiaMath.Frames.Path,Number}",
    "page": "Frames",
    "title": "ModiaMath.Frames.interpolate_r",
    "category": "method",
    "text": "rt = ModiaMath.interpolate_r(path, t)\n\nReturn position r of path::ModiaMath.Path at path parameter t::Number.\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.inverseRotation-Tuple{StaticArrays.SArray{Tuple{3,3},Float64,2,9}}",
    "page": "Frames",
    "title": "ModiaMath.Frames.inverseRotation",
    "category": "method",
    "text": " R_inv = ModiaMath.inverseRotation(R) \n q_inv = ModiaMath.inverseRotation(q)\n\nReturn inverse ModiaMath.RotationMatrixR_inv or inverse ModiaMath.Quaternionq_inv defining the rotation from frame 1 to frame 0  from RotationMatrix R or Quaternion qthat define the rotation from frame 0 to frame 1.\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.planarRotationAngle-Tuple{AbstractArray{T,1} where T,AbstractArray{T,1} where T,AbstractArray{T,1} where T}",
    "page": "Frames",
    "title": "ModiaMath.Frames.planarRotationAngle",
    "category": "method",
    "text": "angle = planarRotationAngle(e, v1, v2; angle_guess = 0.0)\n\nReturn angle of a planar rotation, given the normalized axis of rotation to rotate frame 1 around e into frame 2 (norm(e) == 1 required), and the representations of a vector in frame 1 (v1) and frame 2 (v2). Hereby, it is required that v1 is not parallel to e. The returned angle is in the range -pi <= angle - angle_guess <= pi (from the infinite many solutions, the one is returned that is closest to angle_guess).\n\nExample\n\nimport ModiaMath\nusing Unitful\n\nangle1 = 45u\"°\"\ne      = normalize([1.0, 1.0, 1.0])\nR      = ModiaMath.rot_e(e, angle1)\n\nv1 = [1.0, 2.0, 3.0]\nv2 = ModiaMath.resolve2(R, v1)\n\nangle2 = planarRotationAngle(e, v1, v2)\nisapprox(angle1, angle2)\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.qrot1-Tuple{Number}",
    "page": "Frames",
    "title": "ModiaMath.Frames.qrot1",
    "category": "method",
    "text": "q = ModiaMath.qrot1(angle; q_guess = NullQuaternion)\n\nReturn Quaternion q that rotates with angle angle along the x-axis of frame 1.\n\nFrom the two possible solutions q the one is returned that is closer  to q_guess (note, q and -q define the same rotation).\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.qrot123-Tuple{Number,Number,Number}",
    "page": "Frames",
    "title": "ModiaMath.Frames.qrot123",
    "category": "method",
    "text": "q = ModiaMath.qrot123(angle1, angle2, angle3)\n\nReturn Quaternion q by rotating with angle1 along the x-axis of frame 1, then with angle2 along the y-axis of this frame and then with angle3 along the z-axis of this frame.\n\nFrom the two possible solutions q the one is returned that is closer  to q_guess (note, q and -q define the same rotation).\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.qrot2-Tuple{Number}",
    "page": "Frames",
    "title": "ModiaMath.Frames.qrot2",
    "category": "method",
    "text": "q = ModiaMath.qrot2(angle; q_guess = NullQuaternion)\n\nReturn Quaternion q that rotates with angle angle along the y-axis of frame 1.\n\nFrom the two possible solutions q the one is returned that is closer  to q_guess (note, q and -q define the same rotation).\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.qrot3-Tuple{Number}",
    "page": "Frames",
    "title": "ModiaMath.Frames.qrot3",
    "category": "method",
    "text": "q = ModiaMath.qrot3(angle; q_guess = NullQuaternion)\n\nReturn Quaternion q that rotates with angle angle along the z-axis of frame 1.\n\nFrom the two possible solutions q the one is returned that is closer  to q_guess (note, q and -q define the same rotation).\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.qrot_e-Tuple{StaticArrays.SArray{Tuple{3},Float64,1,3},Number}",
    "page": "Frames",
    "title": "ModiaMath.Frames.qrot_e",
    "category": "method",
    "text": "q = ModiaMath.qrot_e(e, angle; q_guess = NullQuaternion)\n\nReturn Quaternion q that rotates with angle angle along unit axis e. This function assumes that norm(e) == 1.\n\nFrom the two possible solutions q the one is returned that is closer  to q_guess (note, q and -q define the same rotation).\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.qrot_nxy-Tuple{Any,Any}",
    "page": "Frames",
    "title": "ModiaMath.Frames.qrot_nxy",
    "category": "method",
    "text": "q = ModiaMath.qrot_nxy(nx, ny)\n\nIt is assumed that the two input vectors nx and ny are resolved in frame 1 and are directed along the x and y axis of frame 2. The function returns the Quaternion q to rotate from frame 1 to frame 2. \n\nThe function is robust in the sense that it returns always a Quaternion q, even if ny is not orthogonal to nx or if one or both vectors have zero length. This is performed in the following way:  If nx and ny are not orthogonal to each other, first a unit vector ey is  determined that is orthogonal to nx and is lying in the plane spanned by  nx and ny. If nx and ny are parallel or nearly parallel to each other  or ny is a vector with zero or nearly zero length, a vector ey is selected arbitrarily such that ex and ey are orthogonal to each other.  If both nx and ny are vectors with zero or nearly zero length, an arbitrary Quaternion q is returned.\n\nExample\n\nusing Unitful\nimport ModiaMath\n\nq1 = ModiaMath.qrot1(90u\"°\")\nq2 = ModiaMath.qrot_nxy([1  , 0, 0], [0  , 0, 1  ])\nq3 = ModiaMath.qrot_nxy([0.9, 0, 0], [1.1, 0, 1.1])\nisapprox(q1,q2)   # returns true\nisapprox(q1,q3)   # returns true\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.relativeRotation-Tuple{StaticArrays.SArray{Tuple{3,3},Float64,2,9},StaticArrays.SArray{Tuple{3,3},Float64,2,9}}",
    "page": "Frames",
    "title": "ModiaMath.Frames.relativeRotation",
    "category": "method",
    "text": " R_rel = ModiaMath.relativeRotation(R1, R2) \n q_rel = ModiaMath.relativeRotation(q1, q2)\n\nReturn relative ModiaMath.RotationMatrixR_rel or relative ModiaMath.Quaternionq_rel defining the rotation from frame 1 to frame 2  from absolute RotationMatrix R1 or absolute Quaternion q1that define the rotation from frame 0 to frame 1 and the absolute RotationMatrix R2 or the absolute Quaternion q2 that define the rotation from frame 0 to frame 2.\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.resolve1-Tuple{StaticArrays.SArray{Tuple{3,3},Float64,2,9},StaticArrays.SArray{Tuple{3},Float64,1,3}}",
    "page": "Frames",
    "title": "ModiaMath.Frames.resolve1",
    "category": "method",
    "text": "v1 = ModiaMath.resolve1([R|q], v2)\n\nTransform vector v2 (v resolved in frame 2) to vector v1 (v resolved in frame 1) given either ModiaMath.RotationMatrix R or  ModiaMath.Quaternion q (to rotate a frame 1 into a frame 2).\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.resolve2-Tuple{StaticArrays.SArray{Tuple{3,3},Float64,2,9},StaticArrays.SArray{Tuple{3},Float64,1,3}}",
    "page": "Frames",
    "title": "ModiaMath.Frames.resolve2",
    "category": "method",
    "text": "v2 = ModiaMath.resolve2([R|q], v1)\n\nTransform vector v1 (v resolved in frame 1) to vector v2 (v resolved in frame 2) given either ModiaMath.RotationMatrix R or  ModiaMath.Quaternion q (to rotate a frame 1 into a frame 2).\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.rot1-Tuple{Number}",
    "page": "Frames",
    "title": "ModiaMath.Frames.rot1",
    "category": "method",
    "text": "R = ModiaMath.rot1(angle)\n\nReturn RotationMatrix R that rotates with angle angle along the x-axis of frame 1. \n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.rot123-Tuple{Number,Number,Number}",
    "page": "Frames",
    "title": "ModiaMath.Frames.rot123",
    "category": "method",
    "text": "R = ModiaMath.rot123(angle1, angle2, angle3)\n\nReturn RotationMatrix R by rotating with angle1 along the x-axis of frame 1, then with angle2 along the y-axis of this frame and then with angle3 along the z-axis of this frame.\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.rot2-Tuple{Number}",
    "page": "Frames",
    "title": "ModiaMath.Frames.rot2",
    "category": "method",
    "text": "R = ModiaMath.rot2(angle)\n\nReturn RotationMatrix R that rotates with angle angle in [radian] along the y-axis of frame 1. \n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.rot3-Tuple{Number}",
    "page": "Frames",
    "title": "ModiaMath.Frames.rot3",
    "category": "method",
    "text": "R = ModiaMath.rot3(angle)\n\nReturn RotationMatrix R that rotates with angle angle in [radian] along the z-axis of frame 1. \n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.rot_e-Tuple{StaticArrays.SArray{Tuple{3},Float64,1,3},Number}",
    "page": "Frames",
    "title": "ModiaMath.Frames.rot_e",
    "category": "method",
    "text": "R = ModiaMath.rot_e(e, angle)\n\nReturn RotationMatrix that rotates around angle angle along unit axis e. This function assumes that norm(e) == 1.\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.rot_nxy-Tuple{StaticArrays.SArray{Tuple{3},Float64,1,3},StaticArrays.SArray{Tuple{3},Float64,1,3}}",
    "page": "Frames",
    "title": "ModiaMath.Frames.rot_nxy",
    "category": "method",
    "text": "R = ModiaMath.rot_nxy(nx, ny)\n\nIt is assumed that the two input vectors nx and ny are resolved in frame 1 and are directed along the x and y axis of frame 2. The function returns the RotationMatrix R to rotate from frame 1 to frame 2. \n\nThe function is robust in the sense that it returns always a RotationMatrix R, even if ny is not orthogonal to nx or if one or both vectors have zero length. This is performed in the following way:  If nx and ny are not orthogonal to each other, first a unit vector ey is  determined that is orthogonal to nx and is lying in the plane spanned by  nx and ny. If nx and ny are parallel or nearly parallel to each other  or ny is a vector with zero or nearly zero length, a vector ey is selected arbitrarily such that ex and ey are orthogonal to each other.  If both nx and ny are vectors with zero or nearly zero length, an arbitrary rotation matrix is returned.\n\nExample\n\nusing Unitful\nimport ModiaMath\n\nR1 = ModiaMath.rot1(90u\"°\")\nR2 = ModiaMath.rot_nxy([1  , 0, 0], [0  , 0, 1  ])\nR3 = ModiaMath.rot_nxy([0.9, 0, 0], [1.1, 0, 1.1])\nisapprox(R1,R2)   # returns true\nisapprox(R1,R3)   # returns true\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.skew-Tuple{AbstractArray{T,1} where T}",
    "page": "Frames",
    "title": "ModiaMath.Frames.skew",
    "category": "method",
    "text": "M = ModiaMath.skew(e::AbstractVector)\n\nReturn the skew symmetric matrix M::SMatrix{3,3,Float64,9} of vector e (length(e) = 3)\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#ModiaMath.Frames.t_pathEnd-Tuple{ModiaMath.Frames.Path}",
    "page": "Frames",
    "title": "ModiaMath.Frames.t_pathEnd",
    "category": "method",
    "text": "t_end = ModiaMath.t_pathEnd(path::[`ModiaMath.Path`](@ref))\n\nReturn the final path parameter tof the last frame in path (path parameter of first frame = 0.0).\n\n\n\n\n\n"
},

{
    "location": "lib/Frames.html#Frames-1",
    "page": "Frames",
    "title": "Frames",
    "category": "section",
    "text": "Modules = [ModiaMath.Frames]\r\nPrivate = false"
},

{
    "location": "lib/Utilities.html#",
    "page": "Utilities",
    "title": "Utilities",
    "category": "page",
    "text": ""
},

{
    "location": "lib/Utilities.html#ModiaMath.Utilities",
    "page": "Utilities",
    "title": "ModiaMath.Utilities",
    "category": "module",
    "text": "module ModiaMath.Utilities\n\nUtility functions used in ModiaMath.\n\nMain developer\n\nMartin Otter, DLR - Institute of System Dynamics and Control\n\n\n\n\n\n"
},

{
    "location": "lib/Utilities.html#ModiaMath.Utilities.printobj-Tuple{Any,Any}",
    "page": "Utilities",
    "title": "ModiaMath.Utilities.printobj",
    "category": "method",
    "text": "ModiaMath.printobj(name,obj)\n\nPretty print obj as \"<name> = <display(obj)>\". \n\n\n\n\n\n"
},

{
    "location": "lib/Utilities.html#Utilities-1",
    "page": "Utilities",
    "title": "Utilities",
    "category": "section",
    "text": "Modules = [ModiaMath.Utilities]\r\nPrivate = false"
},

{
    "location": "lib/ModiaMath.html#",
    "page": "Constants and Types",
    "title": "Constants and Types",
    "category": "page",
    "text": ""
},

{
    "location": "lib/ModiaMath.html#ModiaMath.AbstractComponentInternal",
    "page": "Constants and Types",
    "title": "ModiaMath.AbstractComponentInternal",
    "category": "type",
    "text": "The internal part of a component (has at least fields \"name\" and \"within\") \n\n\n\n\n\n"
},

{
    "location": "lib/ModiaMath.html#ModiaMath.AbstractComponentWithVariables",
    "page": "Constants and Types",
    "title": "ModiaMath.AbstractComponentWithVariables",
    "category": "type",
    "text": "type ModiaMath.AbstractComponentWithVariables\n\nStruct that contains ModiaMath.AbstractVariables as field or as field in a sub-struct.\n\n\n\n\n\n"
},

{
    "location": "lib/ModiaMath.html#ModiaMath.AbstractRealVariable",
    "page": "Constants and Types",
    "title": "ModiaMath.AbstractRealVariable",
    "category": "type",
    "text": "ModiaMath.AbstractRealVariable <: ModiaMath.AbstractVariable\n\nA real ModiaMath.AbstractVariable (either scalar or array)\n\n\n\n\n\n"
},

{
    "location": "lib/ModiaMath.html#ModiaMath.AbstractSimulationModel",
    "page": "Constants and Types",
    "title": "ModiaMath.AbstractSimulationModel",
    "category": "type",
    "text": "type ModiaMath.AbstractSimulationModel\n\nStruct that is used as simulation model (has field: simulationState)\n\n\n\n\n\n"
},

{
    "location": "lib/ModiaMath.html#ModiaMath.AbstractVariable",
    "page": "Constants and Types",
    "title": "ModiaMath.AbstractVariable",
    "category": "type",
    "text": "type ModiaMath.AbstractVariable <: ModiaMath.AbstractComponentWithVariables\n\nA Variable used as element of the DAE model description and is included in the result (if no residue)\n\n\n\n\n\n"
},

{
    "location": "lib/ModiaMath.html#Constants-and-Types-1",
    "page": "Constants and Types",
    "title": "Constants and Types",
    "category": "section",
    "text": "The following constants are defined in ModiaMath (all these definitions are not exported and therefore need to be prefixed with ModiaMath):const path\nAbsolute path of ModiaMath package directory. This allows for example to run ModiaMath examples as include(\"$(ModiaMath.path)/examples/Simulate_Pendulum.jl\").const Version\nVersion of ModiaMathconst Date\nVersion Date of ModiaMathThe following abstract types are defined and used in ModiaMath:Modules = [ModiaMath]\r\nOrder = [:type]"
},

]}
