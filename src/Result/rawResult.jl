# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.Result (ModiaMath/Result/_module.jl)
#

"""
    mutable struct RawResult
    
Hold result data in a Float64 matrix (first column time, i-th column variable i).
Typical usage:

```
    raw = RawResult(100,3) # Initial storage: 100 time points, 3 variables
    storeRawResult!(raw, [1.0, 2.0, 1.0])
    storeRawResult!(raw, [2.0, 3.0, 5.0])
      ...
    res = getDictResult(res, ["time", "r[1]", "r[2]"]
    
    plot(res["time"], res["r[1]"])
```   
 
If initial storage is not sufficient, it is automatically doubled.
"""
mutable struct RawResult
    nt::Int 
    data::Matrix{Float64}
    names::Vector{String}
      
    function RawResult(nt::Int, names::Vector{String})      
        @assert(nt > 0)
        @assert(length(names) > 0)
        nv = length(names)
        new(0, zeros(nt, nv), names)
    end
end

nResults(res::RawResult) = res.nt

"""
    storeRawResult!(res, v::Vector{Float64})

Store vector `v` in result data structure `res`.
"""
function storeRawResult!(res::RawResult, v::Vector{Float64})
    res.nt += 1
    i = res.nt
    ntMax = size(res.data, 1)

    if i > ntMax
        # Result storage is full, double the current storage
        res.data = [res.data; zeros(ntMax, size(res.data, 2))]
    end
    
    for j = 1:length(v)
        res.data[i,j] = v[j]
    end
    
    return nothing
end

 
function storeRawResult!(res::RawResult, model, t::Float64, x::Vector{Float64}, derx::Vector{Float64}, w::Vector{Float64})
    @assert(length(x) == length(derx))
    res.nt += 1
    i = res.nt
    ntMax = size(res.data, 1)
    
    if i > ntMax
        # Result storage is full, double the current storage
        res.data = [res.data; zeros(ntMax, size(res.data, 2))]
    end

    # Copy data in raw result storage
    res.data[i,1] = t
    jstart = 1
    nx = length(x)
 
    for j = 1:nx
        res.data[i,jstart + j] = x[j]
    end
 
    jstart = jstart + nx
    for j = 1:nx
        res.data[i,jstart + j] = derx[j]
    end
 
    jstart = jstart + nx
    for j = 1:length(w)
        res.data[i,jstart + j] = w[j]
    end
    
    return nothing
end

