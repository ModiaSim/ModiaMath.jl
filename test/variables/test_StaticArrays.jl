"""
    module test_StaticArrays

Evaluate and test properties of StaticArrays
"""
module test_StaticArrays

import ModiaMath

#  Desired:
#    using  StaticArrays
#    import DataFrames
#  
#  In order that these packages need not to be defined in the user environment, they are included via ModiaMath:
using  ModiaMath.StaticArrays
import ModiaMath.DataFrames




function A_add_B!(C, A, B)
    for i = 1:length(A)
        C[i] = A[i] + B[i]
    end
end

function rot(n::AbstractVector, angleRad::Number)
    s  = sin(angleRad)
    c  = cos(angleRad)
    e  = normalize(n)
    es = e * s
    ee = e * e'
    R  = ee + (eye(3) - ee) * c - [ 0.0   -es[3]  es[2];
                                    es[3]  0.0   -es[1];
                                   -es[2]  es[1]  0.0 ]
end


mutable struct ModelVariables
    R1::AbstractMatrix
    R2::MMatrix{3,3,Float64,9}
    R3::SMatrix{3,3,Float64,9}
    R4::SMatrix{3,3,Float64,9}
    v1a::AbstractVector
    v1b::AbstractVector
    v1c::AbstractVector
    v2a::SVector{3,Float64}
    v2b::SVector{3,Float64}
    v2c::SVector{3,Float64}
    v3a::SVector{3,Float64}
    v3b::SVector{3,Float64}
    v3c::SVector{3,Float64}

    function ModelVariables()
        R1  = rot([1,2,-1], 0.5)
        R2  = MMatrix{3,3,Float64,9}(2 * R1)
        R3  = SMatrix{3,3,Float64,9}(3 * R1)
        R4  = SMatrix{3,3,Float64,9}(4 * R1)

        v1a = [1.1, 2.2, 3.3]
        v1b = v1a
        v1c = 3 * copy(v1a)

        v2a = SVector{3,Float64}(v1a)
        v2b = v2a
        v2c = 3 * copy(v2a)

        v3a = SVector{3,Float64}(v1a)
        v3b = v3a
        v3c = 3 * copy(v3a)

        new(R1, R2, R3, R4, v1a, v1b, v1c, v2a, v2b, v2c, v3a, v3b, v3c)
    end
end

const var = ModelVariables()
const nmax = 100000


# Evaluate overhead of for-loop
function test_for_loop(var)
    j = 0
    for i = 1:nmax
        j += 1
    end
    return nothing
end

(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed test_for_loop(var)
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed test_for_loop(var)
println("\ntest_for_loop: elapsedTime = ", elapsedTime, ", bytesAllocated = ", bytesAllocated, "\n")

# Evaluate matrix multiplication and assignment
function mult_StandardArray_assign(var::ModelVariables)::Nothing
    for i = 1:nmax
        var.v1a = var.v1c + i * var.R1 * var.v1c
    end
    return nothing
end


function mult_StandardArray_dotAssign(var::ModelVariables)::Nothing
    for i = 1:nmax
        var.v1a .= var.v1c + i * var.R1 * var.v1c
    end
    return nothing
end

function mult_StandardArray_elementAssign(var::ModelVariables)::Nothing
    for i = 1:nmax
        var.v1a[1:3] = var.v1c + i * var.R1 * var.v1c
    end
    return nothing
end


function mult_StandardArray_mult(var::ModelVariables)::Nothing
    for i = 1:nmax
        A_mul_B!(var.v1a, var.R1, var.v1c)
        A_add_B!(var.v1a, var.v1c, i * var.v1a)
    end
    return nothing
end


function mult_MArray_assign(var::ModelVariables)::Nothing
    for i = 1:nmax
        var.v2a = var.v2c + i * var.R2 * var.v2c
    end
    return nothing
end


function mult_MArray_dotAssign(var::ModelVariables)::Nothing
    for i = 1:nmax
        var.v2a .= var.v2c + i * var.R2 * var.v2c
    end
    return nothing
end


function mult_MArray_elementAssign(var::ModelVariables)::Nothing
    for i = 1:nmax
        var.v2a[1:3] = var.v2c + i * var.R2 * var.v2c
    end
    return nothing
end


function mult_MArray_mult(var::ModelVariables)::Nothing
    for i = 1:nmax
        A_mul_B!(var.v2a, var.R2, var.v2c)
        A_add_B!(var.v2a, var.v2c, i * var.v2a)
    end
    return nothing
end


function mult_SArray_assign(var::ModelVariables)::Nothing
    for i = 1:nmax
        var.v3a = var.v3c + i * var.R3 * var.v3c
    end
    return nothing
end


function mult_SArray_assignWithMMatrix(var::ModelVariables)::Nothing
    for i = 1:nmax
        var.v3a = var.v3c + i * var.R2 * var.v3c
    end
    return nothing
end

#function mult_SArray_dotAssign(var::ModelVariables)
#   for i = 1:nmax
#      var.v3a .= var.R1*var.v3c
#   end
#end


#function mult_SArray_elementAssign(var::ModelVariables)
#   for i = 1:nmax
#      var.v3a[1:3] = var.R1*var.v3c
#   end
#end


table = DataFrames.DataFrame(assignment=String[], relTime=Float64[], elapsedTime=Float64[], bytesAllocated=Int[], link=Bool[])

(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_StandardArray_assign(var)
var.v1b = var.v1a
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_StandardArray_assign(var)
push!(table, ["Vector: va = vc+R*vc", 0.0, elapsedTime, bytesAllocated, var.v1b === var.v1a])

(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_StandardArray_dotAssign(var)
var.v1b = var.v1a
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_StandardArray_dotAssign(var)
push!(table, ["Vector: va .= vc+R*vc", 0.0, elapsedTime, bytesAllocated, var.v1b === var.v1a])

(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_StandardArray_elementAssign(var)
var.v1b = var.v1a
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_StandardArray_elementAssign(var)
push!(table, ["Vector: va[1:3] = vc+R*vc", 0.0, elapsedTime, bytesAllocated, var.v1b === var.v1a])

(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_StandardArray_mult(var)
var.v1b = var.v1a
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_StandardArray_mult(var)
push!(table, ["Vector: A_mult_B!(va,R,vc)", 0.0, elapsedTime, bytesAllocated, var.v1b === var.v1a])


(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_MArray_assign(var)
var.v2b = var.v2a
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_MArray_assign(var)
push!(table, ["SVector: va = vc+R*vc", 0.0, elapsedTime, bytesAllocated, var.v2b === var.v2a])

(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_MArray_dotAssign(var)
var.v2b = var.v2a
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_MArray_dotAssign(var)
push!(table, ["SVector: va .= vc+R*vc", 0.0, elapsedTime, bytesAllocated, var.v2b === var.v2a])

(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_MArray_elementAssign(var)
var.v2b = var.v2a
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_MArray_elementAssign(var)
push!(table, ["SVector: va[1:3] = vc+R*vc", 0.0, elapsedTime, bytesAllocated, var.v2b === var.v2a])

(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_MArray_mult(var)
var.v2b = var.v2a
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_MArray_mult(var)
push!(table, ["SVector: A_mult_B!(va,R,vc)", 0.0, elapsedTime, bytesAllocated, var.v2b === var.v2a])

var.v3b = var.v3a
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_SArray_assignWithMMatrix(var)
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_SArray_assignWithMMatrix(var)
push!(table, ["SVector/MMatrix: va = vc+R*vc", 0.0, elapsedTime, bytesAllocated, var.v3b === var.v3a])

var.v3b = var.v3a
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_SArray_assign(var)
(val, elapsedTime, bytesAllocated, dummy1, dummy2) = @timed mult_SArray_assign(var)
push!(table, ["SVector: va = vc+R*vc", 0.0, elapsedTime, bytesAllocated, var.v3b === var.v3a])


elapsed = table[:,3]
imin = indmin(elapsed)
vmin = elapsed[imin]
for i in 1:length(elapsed)
    table[i,2] = table[i,3] / vmin
end 

println("\nnmax = ", nmax, "\n\ntable = ", table)

end