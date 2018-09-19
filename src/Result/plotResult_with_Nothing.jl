# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.Result (ModiaMath/Result/_module.jl)
#

function plot(result, names::Symbol; heading::AbstractString="", grid::Bool=true, xAxis=:time, figure::Int=1, prefix::AbstractString="", reuse::Bool=false) 
    println("... ModiaMath.plot(..): Call is ignored, since PyPlot not installed in current environment.")
end

function plot(result, names::Tuple; heading::AbstractString="", grid::Bool=true, xAxis=:time, figure::Int=1, prefix::AbstractString="", reuse::Bool=false) 
    println("... ModiaMath.plot(..): Call is ignored, since PyPlot not installed in current environment.")
end

function plot(result, names::AbstractMatrix; heading::AbstractString="", grid::Bool=true, xAxis=:time, figure::Int=1, prefix::AbstractString="", reuse::Bool=false) 
    println("... ModiaMath.plot(..): Call is ignored, since PyPlot not installed in current environment.")
end
