# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module ModiaMath.Result

Organize and plot simulation result data (time series).
The result data of [`ModiaMath.simulate!`](@ref) is returned in one of the
formats supported by this module. The [`ModiaMath.plot`](@ref) function of this
module allows to plot the result data by giving the signal names.
The legends/labels of the plots are automatically constructed by the
signal names and their unit. Example:

```julia
ModiaMath.plot(result, [ (:phi,:r)      (:phi,:phi2,:w);
                         (:w,:w2,:phi2) (:phi,:w)      ],
               heading="Matrix of plots")
```

generates the following plot:

![Matrix-of-Plots](../../resources/images/matrix-of-plots.svg)

# Main developer

[Martin Otter](https://rmc.dlr.de/sr/de/staff/martin.otter/),
[DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)
"""
module Result

export getStringDictResult, SymbolDictResult, StringDictResult, ResultWithVariables
export resultHeading, resultTimeSeries, resultTable, plot
export RawResult, nResults, storeRawResult!
export closeFigure, closeAllFigures
export variablesDependingOnStruct


# using/imports
import ModiaMath
import DataFrames
using  Unitful
using  Requires

# Constants
const headingSize = 10


# include code
include("rawResult.jl")
include("result.jl")
include("plotResult.jl")


# Dummy plot functions if PyPlot is not in the current environment
include("plotResult_with_Nothing.jl")


# If PyPlot is in the current environment, import it in the REPL and install plot functions based on PyPlot
function __init__()
    if !Requires.isprecompiling()
        @eval Main begin
            pyplotAvailable = false
            try
                import PyPlot
                global pyplotAvailable = true
            catch
                println("    PyPlot not available (plot commands will be ignored).\n",
                        "    Try to install PyPlot. See hints here:\n",
                        "    https://github.com/ModiaSim/ModiaMath.jl/wiki/Installing-PyPlot-in-a-robust-way.")
            end

            if pyplotAvailable
                try
                    import PyCall
                catch
                    println("... From ModiaMath: PyCall is added, since needed to set default options for PyPlot")
                    import Pkg
                    Pkg.add("PyCall")
                    import PyCall
                end

                pyplot_rc = PyCall.PyDict(PyPlot.matplotlib."rcParams")
                pyplot_rc["axes.formatter.limits"] = [-3,4]
                pyplot_rc["font.size"]        = 8.0
                pyplot_rc["lines.linewidth"]  = 1.0
                pyplot_rc["grid.linewidth"]   = 0.5
                pyplot_rc["axes.grid"]        = true
                pyplot_rc["axes.titlesize"]   = "medium"
                pyplot_rc["figure.titlesize"] = "medium"

                # Use separate plot windows (no inline plots)
                PyPlot.pygui(true)
            end
        end
    end

    @require PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee" include("plotResult_with_PyPlot.jl")
end


end
