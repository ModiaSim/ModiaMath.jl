# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.Result (ModiaMath/Result/_module.jl)
#


"""
    ModiaMath.plot(result, names; heading="", grid=true, xAxis= :time, 
                   figure=1, prefix="", reuse=false, maxLegend=10)

Plot time series of the result defined by the names keys (Symbol or String).
The keys (and their units, if available in the result) are automatically used as legend.
Units can be either added by using package `Unitful` if result is just a dictionary, or it
can be added by using package `ModiaMath.Result`, where units are defined as elements
of the variable definition. 

# Arguments
Argument `result` maybe one of the following:
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
- `xAxis`: Name of x-axis (Symbol or AbstractString).
- `figure::Int`: Integer identifier of the window in which the diagrams shall be drawn.
- `prefix::AbstractString`: String that is appended in front of every legend label (useful especially if reuse=true)
- `reuse::Bool`: If figure already exists and reuse=false, clear the figure before adding the plot.
- `maxLegend::Int`: If the number of legend entries in one plot command > maxLegend, the legend is suppressed.
  All curves have still their names as labels. The curves can be inspected by their names by clicking in the
  toolbar of the plot on button `Edit axis, curve ..` and then on `Curves`.

# Examples
```julia
import ModiaMath
using Unitful

t = linspace(0.0, 10.0, 100)
result = Dict{AbstractString,Any}(
            "time" => t*u"s", "phi1" => sin.(t)u"rad"  , "phi2" => 0.5*sin.(t),
                              "w1"   => cos.(t)u"rad/s", "w2"   => 0.6*cos.(t))

# 1 signal in one diagram
#   (legend = "phi1 [rad]")
ModiaMath.plot(result, :phi1)   

# 3 signals in one diagram                                 
ModiaMath.plot(result, ("phi1", :phi2, :w1), figure=2)

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

# Append signal of the next simulation run to figure=1
# (legend = "Sim 2: phi1 [rad]")
result[:phi1] = 0.5*result[:phi1]
ModiaMath.plot(result, :phi1, prefix="Sim 2: ", reuse=true)
```

The 5th example above (2 diagrams in form of a vector) give the following plot:

![Figure 5](../../resources/images/plot_figure5.svg)
"""
plot(result, names::Symbol; heading::AbstractString="", grid::Bool=true, xAxis="time", figure::Int=1, prefix::AbstractString="", reuse::Bool=false, maxLegend::Integer=10) =
    plot(result, string(names), heading=heading, grid=grid, xAxis=string(xAxis), figure=figure, prefix=prefix, reuse=reuse, maxLegend=maxLegend)

plot(result, names::AbstractVector; heading::AbstractString="", grid::Bool=true, xAxis="time", figure::Int=1, prefix::AbstractString="", reuse::Bool=false, maxLegend::Integer=10) =
    plot(result, reshape(names, length(names), 1); heading=heading, grid=grid, xAxis=string(xAxis), figure=figure, prefix=prefix, reuse=reuse, maxLegend=maxLegend)
