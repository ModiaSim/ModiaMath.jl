# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control


"""
    module ModiaMath.Variables

Provide variable types for simulation models.

This module provides functions to declare simulation Variables as special structs that
have a `value` and associated attributes. Furthermore, there are functions to copy
the integrator interface variables (`x, derx, residue`) to the appropriate Variables and 
vice versa. Therefore, the modeler does not have to care about, how the values of the variables
are mapped to the integrator interface. Typically, a model is constructed with
macro [`@component`](@ref) using `RealXXX` variable declarations. Example:

```julia
  using ModiaMath
  using StaticArrays

  @component Pendulum(;L=1.0, m=1.0, d=0.1, g=9.81) begin
     phi = RealScalar(start=pi/2, unit="rad"    , fixed=true,               numericType=ModiaMath.XD_EXP)
     w   = RealScalar(start=0.0 , unit="rad/s"  , fixed=true, integral=phi, numericType=ModiaMath.XD_EXP)
     a   = RealScalar(            unit="rad/s^2",             integral=w  , numericType=ModiaMath.DER_XD_EXP) 
     r   = RealSVector{2}(        unit="m"      ,                           numericType=ModiaMath.WC)
  end;
```

The following variable types are currently supported:

| Variable types                       | Description                                    |
|:-------------------------------------|:-----------------------------------------------|
| [`RealScalar`](@ref)                 | Scalar variable with Float64 value             |
| [`RealSVector`](@ref){Size}          | Variable with SVector{Size,Float64} value      |
| [`RealSVector3`](@ref)               | Variable with SVector{3,Float64} value         |
| [`RealVariable`](@ref){VType, EType} | Variable with value/element type `VType/EType` |


All of them have the following attributes:

| Variable attributes                                 | Description                                           |
|:----------------------------------------------------|:------------------------------------------------------|
| `value::VType`                                      | Value of the variable (scale, vector, matrix, ...)    |
| `info::AbstractString`                              | Description text                                      |
| `causality::`[`ModiaMath.Causality`](@ref)          | Causality of variable                                 |
| `variability::`[`ModiaMath.Variability`](@ref)      | Variability of variable                               |
| `start::EType`                                      | Start value of variable                               |
| `fixed::Bool`                                       | `fixed = true`: `value=start` after initialization    |
| `analysis::`[`VariableAnalysisType`](@ref)          | Analysis for which the variable is used               |
| `min::EType`                                        | Minimum value of `value` and of `start`               |
| `max::EType`                                        | Maximum value of `value` and of `start`               |
| `nominal::EType`                                    | Nominal value of `value` (used to improve numerics)   |
| `flow::Bool`                                        | `= true` if variable is a flow variable               |
| `numericyType::`[`ModiaMath.NumericType`](@ref)     | Defines how variable is used by the integrator        |
| `integral`                                          | = the integral variable or `nothing`                  |
| `unit::String`                                      | Unit of `value` and of `start`                        |


The following functions are provided to perform the actual copy-operations and/or to inquire
information about the variables in the model

- [`ModiaMath.copy_start_to_x!`](@ref)
- [`ModiaMath.copy_x_and_derx_to_variables!`](@ref)
- [`ModiaMath.copy_variables_to_residue!`](@ref)
- [`ModiaMath.print_ModelVariables`](@ref)
- [`ModiaMath.get_xTable`](@ref)
- [`ModiaMath.get_copyToVariableTable`](@ref)
- [`ModiaMath.get_copyToResidueTable`](@ref)
- [`ModiaMath.get_copyToResultTable`](@ref)


# Main developer

[Martin Otter](https://rmc.dlr.de/sr/de/staff/martin.otter/), 
[DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)
"""
module Variables

export fullName, instanceName, componentName
export ComponentInternal, isInComponent, isNotInComponent, initComponent!
export @component

export RealVariable, RealScalar, RealSVector, RealSVector3

export ModelVariables
export isArray, isScalar, isUsedInAnalysis, valueLength, indexToString
export get_variableTable, get_xTable, get_copyToVariableTable, get_copyToResidueTable, get_copyToResultTable, print_ModelVariables
export get_variableValueTable, get_residueValueTable
export copy_x_and_derx_to_variables!
export copy_x_to_variables!
export copy_start_to_x!
export copy_variables_to_residue!


export SimulationModel, computeVariables!

# enum values:
export NumericType, Causality, Variability, AnalysisType, VariableAnalysisType
export XD_EXP, XD_IMP, XA, LAMBDA, MUE, DER_XD_EXP, DER_XD_IMP, FD_IMP, FC, WR, WC, TIME, NoNumericType
export Parameter, CalculatedParameter, Input, Output, Local, Independent
export Constant, Fixed, Tunable, Discrete, Continuous
export AnalysisType, KinematicAnalysis, QuasiStaticAnalysis, DynamicAnalysis
export VariableAnalysisType, AllAnalysis, QuasiStaticAndDynamicAnalysis, OnlyDynamicAnalysis, NotUsedInAnalysis

# Enumerations and constants
"""
    @enum NumericType XD_EXP XD_IMP ..

Defines how a variable is used in the integrator interface. The goal is to describe an implicit index-1 DAE:

```math
\\begin{align}
 z &= f_z(x, t) \\\\
 0 &= f_d(\\dot{x}, x, t, z_i > 0) \\\\
 0 &= f_c(x, t, z_i > 0) \\\\
 J &= \\left[ \\frac{\\partial f_d}{\\partial \\dot{x}};  
              \\frac{\\partial f_c}{\\partial x} \\right] \\; \\text{is regular}
\\end{align}
```

The integrator interface to the model is

```julia
getModelResidues(m::AbstractSimulationModel, t::Float64, x::Vector{Float64},  
                 derx::Vector{Float64}, r::Vector{Float64}
```

In the following table it is shown how `NumericType` values of variables are mapped to this
integrator interface:

| @enum value  | Mapping/Description                                                         |
|:-------------|:----------------------------------------------------------------------------|
| `XD_EXP`     | Copied from `x` into variable; der(v) is **computed by model**              |
| `XD_IMP`     | Copied from `x` into variable; der(v) copied from `derx`                    |
| `XA`         | Copied from `x` into variable; `derx` is ignored (algebraic variable)       |
| `DER_XD_EXP` | Computed by model; automatic: `residue = der(v) - derx`                     |
| `DER_XD_IMP` | Copied from `derx` into variable (= der(v); v has type `XD_IMP`             |
| `LAMBDA`     | Copied from `derx` into variable; index-1 algebraic variable                |
| `MUE`        | Copied from `derx`; stabilizing variable to reduced index to index 1        |
| `FD_IMP`     | Is copied from variable to `residue`; part of ``f_d`` equations             |
| `FC`         | Is copied from variable to `residue`; part of ``f_c`` equations             |
| `WR`         | Computed by model; stored in result; used to compute residues               |
| `WC`         | Computed by model; stored in result; only available at communication points |
| `TIME`       | Copied from `t` into variable; independent variable time                    |
"""
@enum NumericType XD_EXP = 1  XD_IMP = 2  XA = 3  LAMBDA = 4  MUE = 5  DER_XD_EXP = 6  DER_XD_IMP = 7 FD_IMP = 8  FC = 9  WR = 10  WC = 11  TIME = 12 NoNumericType = 13


"""
    @enum Causality Parameter CalculatedParameter Input Output Local Independent

Causality of variable (only used when connecting models)
"""
@enum Causality Parameter CalculatedParameter Input Output Local Independent


"""
    @enum Variability Constant Fixed Tunable Discrete Continuous

Currently not used. Will be used in the future to store only the minimum information
in the result (store value only, when it can potentially change).
"""
@enum Variability Constant Fixed Tunable Discrete Continuous


"""
    @enum AnalysisType KinematicAnalysis QuasiStaticAnalysis DynamicAnalysis

Type of analyis that is actually carried out. The `AnalysisType` is set by the user
of the simulation model. Variables are declared in the model with [`VariableAnalysisType`](@ref)`.
Variables with [`VariableAnalysisType`](@ref)` <= AnalysisType` are removed from the analysis and do
not show up in the result. For example, an *acceleration* would be declared as `OnlyDynamicAnalysis`
and then this variable would not show up in the result, if `AnalysisType = KinematicAnalysis` or
`QuasiStaticAnalysis`.
"""
@enum AnalysisType KinematicAnalysis QuasiStaticAnalysis DynamicAnalysis


"""
    @enum VariableAnalysisType AllAnalysis QuasiStaticAndDynamicAnalysis 
                               OnlyDynamicAnalysis NotUsedInAnalysis

Type of analysis that can be carried out with the variable (e.g. a *force* would be defined
as `QuasiStaticAndDynamicAnalysis`, an *acceleration* would be defined as
`OnlyDynamicAnalysis` and a position variable would be defined as
`AllAnalysis`.

Variables with `VariableAnalysisType <= `[`AnalysisType`](@ref) are removed from the analysis and do
not show up in the result. For example, an *acceleration* would be declared as `OnlyDynamicAnalysis`
and then this variable would not show up in the result, if `AnalysisType = KinematicAnalysis` or
`QuasiStaticAnalysis`.
"""
@enum VariableAnalysisType AllAnalysis       QuasiStaticAndDynamicAnalysis OnlyDynamicAnalysis NotUsedInAnalysis

const NumericTypeToVector = [:x, :x, :x, :derx, :derx, :derx, :derx, :residue, :residue, Symbol(""), Symbol(""), Symbol(""), Symbol("")]
const NoNameDefined = :NoNameDefined
const NoUnit = ""

# using/imports
using  StaticArrays
import DataStructures
import DataFrames
import ModiaMath





# include code
include("components.jl")
include("realVariables.jl")
include("modelVariables.jl")
include("simulationModel.jl")

end
