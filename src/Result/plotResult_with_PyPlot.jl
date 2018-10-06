# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.Result (ModiaMath/Result/_module.jl)
#

function prepend!(prefix::AbstractString, ysigLegend::Vector{AbstractString})
   for i in eachindex(ysigLegend)
      ysigLegend[i] = prefix*ysigLegend[i]
   end
   return ysigLegend
end



#--------------------------- Utility plot functions

""" 
    addPlot(result, names, grid, xLabel, xAxis)

Add the time series of one name (if names is one symbol/string) or with
several names (if names is a tuple of symbols/strings) to the current diagram
"""
function addPlot(result, name::AbstractString, grid::Bool, xLabel::Bool, xAxis::AbstractString, prefix::AbstractString, reuse::Bool, maxLegend::Integer)
    (xsig, xsigLegend, ysig, ysigLegend) = resultTimeSeries(result, name, xLabel, xAxis)
    if xsig == nothing; return; end
    if ndims(ysig) == 1
        PyPlot.plot(xsig, ysig, label=prefix*ysigLegend[1])
    else
        for i = 1:size(ysig,2)
            PyPlot.plot(xsig, ysig[:,i], label=prefix*ysigLegend[i])
        end
    end
    PyPlot.grid(grid)
    if length(ysigLegend) <= maxLegend
       PyPlot.legend()
    end
    if xLabel && !reuse
        PyPlot.xlabel(xsigLegend)
    end
end

function addPlot(result, collectionOfNames::Tuple, grid::Bool, xLabel::Bool, xAxis, prefix::AbstractString, reuse::Bool, maxLegend::Integer)
    xsigLegend = ""
    xAxis2 = string(xAxis)
    nLegend = 0

    for name in collectionOfNames
        name2 = string(name)
        (xsig, xsigLegend, ysig, ysigLegend) = resultTimeSeries(result, name2, xLabel, xAxis2)
        if xsig != nothing
            nLegend = nLegend + length(ysigLegend)
            if ndims(ysig) == 1
                PyPlot.plot(xsig, ysig, label=prefix*ysigLegend[1])
            else
                for i = 1:size(ysig,2)
                    PyPlot.plot(xsig, ysig[:,i], label=prefix*ysigLegend[i])
                end            
            end
        end
    end

    PyPlot.grid(grid)
    if nLegend <= maxLegend
       PyPlot.legend()
    end

    if xLabel && !reuse && xsigLegend != nothing 
        PyPlot.xlabel(xsigLegend)
    end
end

addPlot(result, name::Symbol, grid::Bool, xLabel::Bool, xAxis, prefix::AbstractString, reuse::Bool, maxLegend::Integer) =  
        addPlot(result, string(name), grid, xLabel, string(xAxis), prefix, reuse, maxLegend)



#--------------------------- Plot functions
function plot(result, names::AbstractString; heading::AbstractString="", grid::Bool=true, xAxis="time", figure::Int=1, prefix::AbstractString="", reuse::Bool=false, maxLegend::Integer=10) 
    PyPlot.figure(figure) 
    if !reuse
       PyPlot.clf()
    end
    addPlot(result, names, grid, true, string(xAxis), prefix, reuse, maxLegend)
    heading2 = getHeading(result, heading)
    
    if heading2 != "" && !reuse
        PyPlot.title(heading2, size=headingSize)    # Python 2.7: fontsize; Python 3.x: size
      #PyPlot.suptitle(heading, fontsize=12)
      #PyPlot.tight_layout()
      #PyPlot.subplots_adjust(top=0.88)
    end
end

function plot(result, names::Tuple; heading::AbstractString="", grid::Bool=true, xAxis="time", figure::Int=1, prefix::AbstractString="", reuse::Bool=false, maxLegend::Integer=10) 
    PyPlot.figure(figure)
    if !reuse
       PyPlot.clf()
    end
    addPlot(result, names, grid, true, string(xAxis), prefix, reuse, maxLegend)

    heading2 = getHeading(result, heading)
    
    if heading2 != "" && !reuse
        PyPlot.title(heading2, size=headingSize)    # Python 2.7: fontsize; Python 3.x: size
    end
end

function plot(result, names::AbstractMatrix; heading::AbstractString="", grid::Bool=true, xAxis="time", figure::Int=1, prefix::AbstractString="", reuse::Bool=false, maxLegend::Integer=10) 
    xAxis2 = string(xAxis)
    PyPlot.figure(figure)
    if !reuse
       PyPlot.clf()
    end
    heading2 = getHeading(result, heading)
    (nrow, ncol) = size(names)    

    # Add signals
    k = 1
    for i = 1:nrow
        xLabel = i == nrow
        for j = 1:ncol
            PyPlot.subplot(nrow, ncol, k)
            addPlot(result, names[i,j], grid, xLabel, xAxis2, prefix, reuse, maxLegend)
            k = k + 1
            if ncol == 1 && i == 1 && heading2 != "" && !reuse
                PyPlot.title(heading2, size=headingSize)
            end
        end
    end

    # Add overall heading in case of a matrix of diagrams (ncol > 1)
    if ncol > 1 && heading2 != "" && !reuse
        PyPlot.suptitle(heading2, size=headingSize)
        #   PyPlot.tight_layout()
        #   PyPlot.subplots_adjust(top=0.88)
    end
end
