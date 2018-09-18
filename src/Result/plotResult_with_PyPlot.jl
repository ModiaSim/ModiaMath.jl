# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.Result (ModiaMath/Result/_module.jl)
#


#--------------------------- Utility plot functions
""" 
    addPlot(result, names, grid, xLabel, xAxis)

Add the time series of one name (if names is one symbol/string) or with
several names (if names is a tuple of symbols/strings) to the current diagram
"""
function addPlot(result, name::Symbol, grid::Bool, xLabel::Bool, xAxis)
    (xsig, xsigLegend, ysig, ysigLegend) = resultTimeSeries(result, name, xLabel, xAxis)
    if xsig == nothing; return; end
    PyPlot.plot(xsig, ysig)
    PyPlot.grid(grid)
    PyPlot.legend(ysigLegend)
    if xLabel
        PyPlot.xlabel(xsigLegend)
    end
end

function addPlot(result, collectionOfNames::Tuple, grid::Bool, xLabel::Bool, xAxis)
    legend = String[]
    xsigLegend = ""

    for name in collectionOfNames
        name2 = Symbol(name)
        (xsig, xsigLegend, ysig, ysigLegend) = resultTimeSeries(result, name2, xLabel, xAxis)
        if xsig != nothing
            PyPlot.plot(xsig, ysig)
            append!(legend, ysigLegend)
        end
    end

    PyPlot.grid(grid)
    PyPlot.legend(legend)

    if xLabel && xsigLegend != nothing
        PyPlot.xlabel(xsigLegend)
    end
end

addPlot(result, name::String, grid::Bool, xLabel::Bool, xAxis) =  addPlot(result, Symbol(name), grid, xLabel, xAxis)



#--------------------------- Plot functions
function plot(result, names::Symbol; heading::AbstractString="", grid::Bool=true, xAxis=:time, figure::Int=1) 
    PyPlot.figure(figure)
    PyPlot.clf()
    addPlot(result, names, grid, true, xAxis)
    heading2 = getHeading(result, heading)
    
    if heading2 != ""
        PyPlot.title(heading2, size=headingSize)    # Python 2.7: fontsize; Python 3.x: size
      #PyPlot.suptitle(heading, fontsize=12)
      #PyPlot.tight_layout()
      #PyPlot.subplots_adjust(top=0.88)
    end
end

function plot(result, names::Tuple; heading::AbstractString="", grid::Bool=true, xAxis=:time, figure::Int=1) 
    PyPlot.figure(figure)
    PyPlot.clf()
    addPlot(result, names, grid, true, xAxis)

    heading2 = getHeading(result, heading)
    
    if heading2 != ""
        PyPlot.title(heading2, size=headingSize)    # Python 2.7: fontsize; Python 3.x: size
    end
end

function plot(result, names::AbstractMatrix; heading::AbstractString="", grid::Bool=true, xAxis=:time, figure::Int=1) 
    PyPlot.figure(figure)
    PyPlot.clf()
    heading2 = getHeading(result, heading)
    (nrow, ncol) = size(names)    

    # Add signals
    k = 1
    for i = 1:nrow
        xLabel = i == nrow
        for j = 1:ncol
            PyPlot.subplot(nrow, ncol, k)
            addPlot(result, names[i,j], grid, xLabel, xAxis)
            k = k + 1
            if ncol == 1 && i == 1 && heading2 != ""
                PyPlot.title(heading2, size=headingSize)
            end
        end
    end

    # Add overall heading in case of a matrix of diagrams (ncol > 1)
    if ncol > 1 && heading2 != ""
        PyPlot.suptitle(heading2, size=headingSize)
        #   PyPlot.tight_layout()
        #   PyPlot.subplots_adjust(top=0.88)
    end
end
