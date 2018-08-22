# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.Result (ModiaMath/Result/_module.jl)
#


const headingSize = 10

""" 
    addPlot(result, names, grid, xLabel, xAxis)

Add the time series of one name (if names is one symbol/string) or with
several names (if names is a tuple of symbols/strings) to the current diagram
"""
addPlot(result, name::String, grid::Bool, xLabel::Bool, xAxis) =  addPlot(result, Symbol(name), grid, xLabel, xAxis)



"""
    ModiaMath.plot(result, names; heading="", grid=true, xAxis= :time, figure=1)

Plot time series of the result defined by the names keys (Symbol or String).
The keys (and their units, if available in the result) are automatically used as legend.
Units can be either added by using package `Unitful` if result is just a dictionary, or it
can be added by using package `ModiaMath.Result`, where units are defined as elements
of the variable definition. 

# Arguments
Argument `result` maybe one of the following:
- A dictionary `Dict{Symbol,AbstractVector}`
- A dictionary `Dict{String,AbstractVector}`
- A dictionary `Dict{AbstractString,Any}`
- An instance of struct [`ModiaMath.Result`](@ref)
- An object for which function [`ModiaMath.resultTimeSeries`](@ref) is defined.

Argument `names` defines the diagrams to be drawn and the time series to be included in the respective diagram: 
- If names is a **Symbol** or **String**, generate one diagram with one time series.
- If names is a **Tuple** of Symbols/Strings, generate one diagram with the time series of the given keys
- If names is a **Vector** or **Matrix** of Symbols/Strings/Tuples, generate a vector or matrix of diagrams.

Remaining arguments:
- `heading::AbstractString`: Optional heading above the diagram.
- `grid::Bool`: Optional grid.
- `xAxis`: Name of x-axis (Symbol or String).
- `figure::Int`: Integer identifier of the window in which the diagrams shall be drawn.

# Examples
```julia
import ModiaMath
using Unitful

t = linspace(0.0, 10.0, 100)
result = Dict{Symbol,Vector{Float64}}(
            :time=>t*u"s", :phi1=>sin.(t)u"rad", :phi2=>0.5*sin.(t),
                           :w1  =>cos.(t)u"rad/s", :w2  => 0.6*cos.(t))

# 1 signal in one diagram
ModiaMath.plot(result, :phi1)

# 3 signals in one diagram                                 
ModiaMath.plot(result, (:phi1, :phi2, :w1), figure=2)

# 3 diagrams in form of a vector (every diagram has one signal)                 
ModiaMath.plot(result, [:phi1, :phi2, :w1], figure=3)     

# 4 diagrams in form of a matrix (every diagram has one signal)          
ModiaMath.plot(result, ["phi1" "phi2";
                        "w1"   "w2"   ], figure=4)     

# 2 diagrams in form of a vector           
ModiaMath.plot(result, [ (:phi1,:phi2), (:w1) ], figure=5)           

# 4 diagrams in form of a matrix
ModiaMath.plot(result, [ (:phi1,)           (:phi2,:w1);
                         (:phi1,:phi2,:w1)  (:w2,)     ],figure=6)  

# Plot w1=f(phi1) in one diagram 
ModiaMath.plot(result, :w1, xAxis=:phi1, figure=7)                   
```

The 5th example above (2 diagrams in form of a vector) give the following plot:

![Figure 5](../../resources/images/plot_figure5.svg)
"""
plot(result, names::String; heading="", grid=true, xAxis=:time, figure=1) =
    plot(result, Symbol(names), heading=heading, grid=grid, xAxis=xAxis, figure=figure)

plot(result, names::AbstractVector; heading::AbstractString="", grid::Bool=true, xAxis=:time, figure::Int=1) =
    plot(result, reshape(names, length(names), 1); heading=heading, grid=grid, xAxis=xAxis, figure=figure)
