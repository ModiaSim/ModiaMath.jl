 # License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.Result (ModiaMath/Result/_module.jl)
#

#---------------------------- Acessing the time series results (various options)
appendUnit(name, unit) = unit == "" ? string(name) : string(name, " [", unit, "]")

function getStringDictResult(model,res::RawResult)
   nt     = res.nt
   data   = res.data
   result = Dict{String,AbstractVector}()
   for (i,name) in enumerate(res.names)
      result[name] = @view data[1:nt,i]
   end
   return result, nt
end

const SymbolDictResult = Dict{Symbol,AbstractVector}
resultHeading(result::SymbolDictResult) = ""

"""
    (xsig, xsigLegend, ysig, ysigLegend) = 
           ModiaMath.resultTimeSeries(result, name, xLabel::Bool, xAxis)

For a desired `result` data structure, this function has to be provided
to return the x-vector (`xsig`), the y-vector (`ysig`) and the legend
of the x-vector (`xsigLegend`), and of the y-vector (`ysigLegend`) as Strings, given 
the key of the y-vector (`name`) and the key of the x-vector (`xAxis`).
If `xLabel=false` the legend of the x-vector should be an empty string (`""`).
"""
function resultTimeSeries(result::SymbolDictResult, name, xLabel::Bool, xAxis)
    xsig = get(result, Symbol(xAxis), nothing)
    if xsig == nothing
       warn("\nModiaMath.plot: argument xAxis (= ", xAxis, ") is not correct or does not identify a signal in the result.")
       return (nothing, nothing, nothing, nothing)
    end
    xsigLegend = xLabel ? appendUnit( string(xAxis), string(unit(xsig[1])) ) : ""
    xsig       = ustrip.( xsig )

    ysig = get(result, Symbol(name), nothing)
    if ysig == nothing
       warn("\nModiaMath.plot: argument name (= ", name, ") is not correct or does not identify a signal in the result.")
       return (nothing, nothing, nothing, nothing)
    end 
    ysigLegend = [appendUnit( string(name), string(unit(ysig[1])) )]
    ysig       = ustrip.( ysig )

    return(xsig, xsigLegend, ysig, ysigLegend)
end


const StringDictResult = Dict{String,AbstractVector}
resultHeading(result::StringDictResult) = ""
function resultTimeSeries(result::StringDictResult, name, xLabel::Bool, xAxis)
    xsig = get(result, string(xAxis), nothing)
    if xsig == nothing
       warn("\nModiaMath.plot: argument xAxis (= ", xAxis, ") is not correct or does not identify a signal in the result.")
       return (nothing, nothing, nothing, nothing)
    end
    xsigLegend = xLabel ? appendUnit( string(xAxis), string(unit(xsig[1])) ) : ""
    xsig       = ustrip.( xsig )

    ysig = get(result, string(name), nothing)
    if ysig == nothing
       warn("\nModiaMath.plot: argument name (= ", name, ") is not correct or does not identify a signal in the result.")
       return (nothing, nothing, nothing, nothing)
    end 
    ysigLegend = [appendUnit( string(name), string(unit(ysig[1])) )]
    ysig       = ustrip.( ysig )

    return(xsig, xsigLegend, ysig, ysigLegend)
end


const StringDictAnyResult = Dict{AbstractString,Any}
resultHeading(result::StringDictAnyResult) = ""
function resultTimeSeries(result::StringDictAnyResult, name, xLabel::Bool, xAxis)
    xsig = get(result, string(xAxis), nothing)
    if xsig == nothing
       warn("\nModiaMath.plot: argument xAxis (= ", xAxis, ") is not correct or does not identify a signal in the result.")
       return (nothing, nothing, nothing, nothing)
    end
    xsigLegend = xLabel ? appendUnit( string(xAxis), string(unit(xsig[1])) ) : ""
    xsig       = ustrip.( xsig )

    ysig = get(result, string(name), nothing)
    if ysig == nothing
       warn("\nModiaMath.plot: argument name (= ", name, ") is not correct or does not identify a signal in the result.")
       return (nothing, nothing, nothing, nothing)
    end 
    ysigLegend = [appendUnit( string(name), string(unit(ysig[1])) )]
    ysig       = ustrip.( ysig )

    return(xsig, xsigLegend, ysig, ysigLegend)
end


"""
    mutable struct ResultWithVariables

Struct that is generated as a result of ModiaMath.simulate!(...),
if the model is a Modia3D.AbstractComponentWithVariables struct.
"""
mutable struct ResultWithVariables
   "Dictionary of time series (elements can be Float64/Int64/Bool Vectors or Matrices)"
   series::Dict{Symbol,Union{AbstractVector,AbstractMatrix}}  

   "Dictionary of variables (elements are <: ModiaMath.AbstractVariable with required fields value and unit)"
   var::Dict{Symbol,Any}   

   "String used as optional heading of a plot window"                                   
   resultHeading::String

   ResultWithVariables(series, var, resultHeading="") = new(series, var, resultHeading)
end

resultHeading(result::ResultWithVariables) = result.resultHeading


function getSignal(seriesDict, name)
   nameAsString = string(name)
   if haskey(seriesDict, name)
      sig     = seriesDict[name]
      keyName = name
   else
      if nameAsString[end] == ']'
         indexRange = rsearch(nameAsString, "[")
         i = indexRange[1]
         @assert(i >= 2)
         keyName = Symbol(nameAsString[1:i-1])
 
         if haskey(seriesDict, keyName)
            sig2 = seriesDict[keyName]

            # Build AST and evaluate: sig = sig2[:,...]
            ex1 = parse(nameAsString[i:end])
            ex2 = :( $sig2[:] )
            append!(ex2.args, ex1.args)
            sig = @eval $ex2
         else
            warn("\nModiaMath.plot: argument name (= ", name, ") is not correct.")
            sig = nothing
         end  
      else
         warn("\nModiaMath.plot: argument name (= ", name, ") is not correct or does not identify a signal in the result.")     
         sig = nothing
         keyName = name
      end
   end
   return (sig, keyName, nameAsString)
end


function resultTimeSeries(result::ResultWithVariables, name, xLabel::Bool, xAxis)
   seriesDict = result.series
   (ysig,ykeyName,yNameAsString) = getSignal(seriesDict, Symbol(name))
   if ysig == nothing
      return (nothing, nothing, nothing, nothing) 
   end
   yvar = result.var[ykeyName]
   if ndims(ysig) == 1
      # ysig is a vector
      ysigLegend = [appendUnit(yNameAsString, yvar.unit)]
   else
      # sig has more as one dimension
      ysigLegend = Array{String}(length(yvar.value))
      for i in eachindex(yvar.value)
         ysigLegend[i] = appendUnit(ModiaMath.indexToString(ykeyName, yvar.value, i), yvar.unit)
      end
   end
   (xsig,xkeyName,xNameAsString) = getSignal(seriesDict, Symbol(xAxis))
   xsigLegend = xLabel ? appendUnit(xNameAsString, result.var[xkeyName].unit) : ""
   return (xsig, xsigLegend, ysig, ysigLegend)
end


getHeading(result, heading) = heading!="" ? heading : resultHeading(result)

