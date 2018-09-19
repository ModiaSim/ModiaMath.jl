# Plans for version 1.0

ModiaMath is not yet ready and should not be used for production simulation runs.
The following features are planned to be implemented before ModiaMath version 1.0:

- **General**
  - Support of all features needed by Modia 1.0 and Modia3D 1.0.
  - Improved documentation.

- **SimulationEngine**
  - Sparse Jacobian for IDA (mainly available, but not yet transferred to restructured simulation engine).
  - Sparse Jacobian for nonlinear solver during initialization.
  - User-supplied dense Jacobian for IDA (pre-requisite for next step).
  - Support for mue-stabilization technique for DAEs with index > 1
    (using Jacobian information).
  - Change basic engine from IDA to package DifferentialEquations (and access IDA via this interface).
  - If model is an ODE, support ODE integrators from DifferentialEquations package.
  - New inquiry function isCompletedIntegratorStep(..) to inquire when the integrator step 
    is completed. This function shall have the same functionality as fmi2CompletedIntegratorStep and
    requires to restructure the integrator while-loop since the integrator must return after every
    completed step.
  - Linearization of the DAE after initialization and at the actual integration time instant in order
    that linear analysis and synthesis methods can be applied.

- **NonlinearEquations**
  - Transfering the FORTRAN code hybrd.f from Minpack.
    to a native Julia implementation. It is expected that this solver
    is more robust as Sundials KINSOL (implemented with adaptations, but not yet fully tested).
  - Using Minpack during initialization.

- **Result**
  - Integer and Boolean result variables.
  - Several time vectors for different result data.
  - Parameters (time vector has length 2).
  - Clocked variables (time vector has only values at clock ticks).
  - Automatic plot layout for clocked variables.
  - Support of another plot package (besides PyPlot). Probably best to 
    support Plots.
  - Use named tuples to define line style elements, e.g.\
    `plot(result, ((name=:phi1, color=:blue), (name=:phi2, color=:red)) )`

- **Frames**
  - Support splines for the interpolation of Frames
    (currently only linear interpolation supported).

