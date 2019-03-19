# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.Variables (ModiaMath/Variables/_module.jl)
#

# Definitions of Real variables

#=
The dynamic system is described as index-1 DAE:

    wr := wr(der(x),x,t)      # Additional variables, that are explicitly computed but are not reported to the integrator (if used in fc(...), wr is not allowed to depend on der(x)
    0   = fd(der(x),x,t,wr)   
    0   = fc(x,t,wr) 
    wc := wc(der(x),x,t,wr)   # Additional variables, that need to be only computed at communication points.

such that (assuming wr is inlined in fd/fc)
  
   [der(fd, der(x));
    der(fc, x)     ] is regular   

The vectors are partitioned as:

   x      = [xd_exp;        # XD_EXP: variables appearing differentiated in the model and the model explicitly computes der(xd_exp), so der(xd_exp) is an output of the model (used in model)
             xd_imp;        # XD_IMP: variables appearing differentiated in the model and the model does not compute der(xd_imp), so der(xd_imp) is an input to the model (used in model)
             xa;            # XA    : variables not appearing differentiated in the model, so algebraic variables (used in model)
             lambda_int;    #         integral of index-1 algebraic variables (not used in model; used only by integrator)
             mue_int   ]    #         integral of stabilizing variables der(mue_int) (not used in model; used only by integrator)

   der(x) = [der(xd_exp);   # DER_XD_EXP: derivative of xd_exp; is computed explicitly in the model (used in model).
             der(xd_imp);   # DER_XD_IMP: derivative of xd_imp; is provided by integrator (used in model)
             der(xa);       #             derivative of xa (not used in model; used only in integrator)
             lambda;        # LAMBDA    : index-1 algebraic variable; is provided by integrator (used in model)
             mue        ]   # MUE       : stabilizing variable to reduce higher index to index 1; is provided by integrator (used in model)

   residue = [fd_exp;       #           : residue of explicitly solved fd equations; computed by code below (e.g. if x[1] = phi and x[2] = v (= der(phi)), then residue[1] = derx[1] - v)
              fd_imp;       # FD_IMP    : residue of implicitly solved fd equations (is provided by model)
              fc]           # FC        : residue of fc equations (is provided by model)

If x consists only of xd_exp, then the DAE is an ODE and an ODE integration method can be used.

All variables x and der(x) are subtypes of type Real (typically Float64).
These variables are below defined as subtypes of type ContinuousTimeVariable

Types WR are auxiliary variables that are explicitly computed as functions of t,x,der(x) and are used to compute residues
(so these variables are computed at every call of the model)

Types WC are additional variables that are computed only at communication points.

Type TIME is just used to mark which variable is the independent time variable.
=#

isNothing(T) = typeof(T) == Nothing
hasValue(T)  = typeof(T) != Nothing


#------------------------------ Generic Real Variable ----------------------------------------
"""
    v = RealVariable{ValueType,ElementType}(...)

Generate a new variable with `v.value::ValueType`, `v.nominal::ElementType`.
The argument list is described in module [`Variables`](@ref).
"""
mutable struct RealVariable{ValueType,ElementType} <: ModiaMath.AbstractRealVariable
    # Generic Attributes of every AbstractVariable (partially based on FMI 2.0)
    _internal::ComponentInternal
    value::ValueType                 # Actual value of variable (can be a scalar or an array)
    info::AbstractString
    causality::Causality
    variability::Variability
    start::ValueType                 # Initial value of variable
    fixed::Bool                      # = true, start is fixed during initialization; = false, start is a guess value (can be changed during initialization)
    analysis::VariableAnalysisType   # Analysis for which the Variable is used

    # Attributes specific to Real variables
    min::ElementType
    max::ElementType
    nominal::ElementType                                              # nominal value; is used to compute absolute tolerances and might be used for scaling
    flow::Bool
    numericType::NumericType                                          # how the variable is used in the equations
    integral::Union{Nothing,RealVariable{ValueType,ElementType}}      # if present, integral is the variable that represents the integral of the actual variable (so variable = d(integral)/dt     derivative::Union{Nothing, RealVariable{ValueType, ElementType}}  # if present, derivative is the variable that represents the derivative of the actual variable (so derivative = d(variable)/dt)
    derivative::Union{Nothing,RealVariable{ValueType,ElementType}}    # if present, derivative is the variable that represents the derivative of the actual variable (so derivative = d(variable)/dt)
    unit::String                                                      # unit of the variable (temporal solution until package Unitful is supported in Julia v0.7)

    # How the variable is stored in vectors
    ivar::Int        # if value is stored in x/derx/residue vector: x/derx/residue[ ivar ] = first value of variable.value
    iresult::Int     # if value is stored in result vector: result[ iresult ] = first value of variable.value

    function RealVariable{ValueType,ElementType}(name, within, info, start, unit, fixed, min, max, nominal, flow, causality, variability, 
      numericType, integral, analysis) where {ValueType, ElementType}

        variable = new(ComponentInternal(name, within), deepcopy(start), info, causality, variability, 
                     start, fixed, analysis, min, max, nominal, flow, numericType, integral, nothing, unit, 0, 0) 

        if typeof(within) != Nothing
            setfield!(within, Symbol(name), variable)
        end

        if hasValue(integral)
            integral.derivative = variable      
        end

        return variable      
    end
end

RealVariable{ValueType,ElementType}(name=NoNameDefined,
    within=nothing;
    info="",
    start=nothing,
    unit=NoUnit,
    fixed=false,
    min=-ElementType(Inf),
    max=ElementType(Inf),
    nominal=NaN, 
    flow=false,
    causality=Local, 
    variability=Continuous, 
    integral=nothing, 
    numericType=WR,
    analysis=ModiaMath.AllAnalysis) where {ValueType, ElementType} =
        RealVariable{ValueType,ElementType}(name, within, info, start, unit, fixed, min, max, nominal, false, causality, variability, numericType, integral, analysis)

function Base.show(io::IO, r::RealVariable{ValueType,ElementType}) where {ValueType, ElementType}
    print(r.value, " ", r.unit)
    if r.info != ""
      # print(" # ", r.info)
   end
end



#------------------------------ Concrete Real Variable types ----------------------------------------

"""
    v = RealScalar(..)

Generate a variable of type [`RealVariable`](@ref)`{Float64,Float64}` to define scalar, Float64, real variables.
The argument list is described in module [`Variables`](@ref).
"""
const RealScalar = RealVariable{Float64,Float64}


"""
    v = RealSVector{Size}(..)

Generate a variable of type [`RealVariable`](@ref)`{SVector{Size,Float64},Float64}` where the values have
type StaticArrays.SVector{Size,Float64}.
The argument list is described in module [`Variables`](@ref).
"""
const RealSVector{Size} = RealVariable{SVector{Size,Float64},Float64}


"""
    v = RealSVector3(..)

Generate a variable of type [`RealVariable`](@ref)`{SVector{3,Float64},Float64}` where the values have
type StaticArrays.SVector{3,Float64}.
The argument list is described in module [`Variables`](@ref).
"""
const RealSVector3 = RealSVector{3}


RealScalar(name=NoNameDefined,
    within=nothing;
    info="",
    start=0.0,
    unit=NoUnit,
    fixed=false,
    min=-Inf,
    max=Inf,
    nominal=NaN, 
    flow=false,
    causality=Local, 
    variability=Continuous, 
    integral=nothing, 
    numericType=WR,
    analysis=ModiaMath.AllAnalysis) =
         RealScalar(name, within, info, start, unit, fixed, min, max, nominal, flow, causality, variability, numericType, integral, analysis)

RealSVector{Size}(name=NoNameDefined,
    within=nothing;
    info="",
    start=zeros(Size),
    unit=NoUnit,
    fixed=false,
    min=-Inf,
    max=Inf,
    nominal=NaN, 
    flow=false,
    causality=Local, 
    variability=Continuous, 
    integral=nothing, 
    numericType=WR,
    analysis=ModiaMath.AllAnalysis) where {Size} =
        RealSVector{Size}(name, within, info, start, unit, fixed, min, max, nominal, flow, causality, variability, numericType, integral, analysis)

