# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module
#   ModiaMath.Frames (ModiaMath/Frames/_module.jl)
#


"""
    const ModiaMath.Quaternion = SVector{4,Float64}

Describes the rotation from a frame 1 into a frame 2 with a quaternion vector.
If `e` is the (normalized) axis of rotation to rotate frame 1 into frame 2
(either resolved in frame 1 or frame 2) and `angle`
is the rotation angle for this rotation then the quaternion vector
`q::ModiaMath.Quaternions` is defined as:

```julia
q = [e*sin(angle/2),
       cos(angle/2]
```
"""
const Quaternion = SVector{4,Float64}


"""
    const ModiaMath.NullQuaternion = Quaternion(0,0,0,1)

Constant Quaternion vector of a null rotation (= no rotation from frame 1 to frame 2)
"""
NullQuaternion = Quaternion(0.0, 0.0, 0.0, 1.0)


"""
    ModiaMath.assertQuaternion(q::AbstractVector)

Assert that vector `q` has the properties of a `Quaternion` vector
(has 4 elements, `norm(q) = 1`)
"""
function assertQuaternion(q::AbstractVector)
    @assert(length(q) == 4)
    @assert(abs(norm(q) - 1.0) <= 1e-10)
end


"""
    R = ModiaMath.from_q(q::ModiaMath.Quaternion)

Return RotationMatrix `R` from Quaternion `q`.
"""
from_q(q::Quaternion) = RotationMatrix(2.0 * (q[1] * q[1] + q[4] * q[4]) - 1.0, 
                                       2.0 * (q[2] * q[1] - q[3] * q[4]),
                                       2.0 * (q[3] * q[1] + q[2] * q[4]),
                                       2.0 * (q[1] * q[2] + q[3] * q[4]),
                                       2.0 * (q[2] * q[2] + q[4] * q[4]) - 1.0,
                                       2.0 * (q[3] * q[2] - q[1] * q[4]),
                                       2.0 * (q[1] * q[3] - q[2] * q[4]), 
                                       2.0 * (q[2] * q[3] + q[1] * q[4]), 
                                       2.0 * (q[3] * q[3] + q[4] * q[4]) - 1.0)



const p4limit = 0.1
const c4limit = 4.0 * p4limit * p4limit

"""
    q = ModiaMath.from_R(R::ModiaMath.RotationMatrix;
                         q_guess = NullQuaternion)

Return `Quaternion q` from `RotationMatrix R`.

From the two possible solutions `q` the one is returned that is closer 
to `q_guess` (note, `q` and `-q` define the same rotation).
"""
function from_R(R::RotationMatrix;
                q_guess::Quaternion=NullQuaternion)::Quaternion
   #= 
   Below it is guaranteed that c1>=0, c2>=0, c3>=0, c4>=0 and
   that not all of them can be zero at the same time
   (e.g., if 3 of them are zero, the 4th variable is 1).
   Since the sqrt(..) has to be performed on one of these variables,
   it is applied on a variable which is far enough from zero.
   This guarantees that the sqrt(..) is never taken near zero
   and therefore the derivative of sqrt(..) can never be infinity.
   There is an ambiguity for quaternions, since q and -q
   lead to the same RotationMatrix. This ambiguity
   is resolved here by selecting the q that is closer to
   the input argument q_guess.   
   =#
    c1::Float64 = 1.0 + R[1,1] - R[2,2] - R[3,3]
    c2::Float64 = 1.0 + R[2,2] - R[1,1] - R[3,3]
    c3::Float64 = 1.0 + R[3,3] - R[1,1] - R[2,2]
    c4::Float64 = 1.0 + R[1,1] + R[2,2] + R[3,3]
    paux::Float64  = 0.0
    paux4::Float64 = 0.0
 
    if c4 > c4limit || (c4 > c1 && c4 > c2 && c4 > c3) 
        paux  = sqrt(c4) / 2
        paux4 = 4 * paux
        q = Quaternion((R[2,3] - R[3,2]) / paux4,
                       (R[3,1] - R[1,3]) / paux4,
                       (R[1,2] - R[2,1]) / paux4,
                       paux)
   
    elseif c1 > c2 && c1 > c3 && c1 > c4 
        paux  = sqrt(c1) / 2
        paux4 = 4 * paux
        q = Quaternion(paux,
                    (R[1,2] + R[2,1]) / paux4,
                    (R[1,3] + R[3,1]) / paux4,
                    (R[2,3] - R[3,2]) / paux4)
   
    elseif c2 > c1 && c2 > c3 && c2 > c4 
        paux  = sqrt(c2) / 2
        paux4 = 4 * paux
        q = Quaternion((R[1,2] + R[2,1]) / paux4,
                      paux,
                     (R[2,3] + R[3,2]) / paux4,
                     (R[3,1] - R[1,3]) / paux4)
   
    else
        paux  = sqrt(c3) / 2
        paux4 = 4 * paux
        q = Quaternion((R[1,3] + R[3,1]) / paux4,
                     (R[2,3] + R[3,2]) / paux4,
                      paux,
                     (R[1,2] - R[2,1]) / paux4)
    end
   
    return dot(q, q_guess) >= 0 ? q : -q
end

from_R(R::AbstractMatrix; q_guess::AbstractVector=NullQuaternion)::Quaternion = 
    from_R(RotationMatrix(R); q_guess=Quaternion(q_guess))
    


"""
    q = ModiaMath.qrot1(angle; q_guess = NullQuaternion)

Return Quaternion `q` that rotates with angle `angle` along the x-axis of frame 1.

From the two possible solutions `q` the one is returned that is closer 
to `q_guess` (note, `q` and `-q` define the same rotation).
"""
@inline function qrot1(angle::Number; q_guess::Quaternion=NullQuaternion)::Quaternion
    q = Quaternion(sin(angle / 2), 0.0, 0.0, cos(angle / 2))

    return dot(q, q_guess) >= 0 ? q : -q
end


"""
    q = ModiaMath.qrot2(angle; q_guess = NullQuaternion)

Return Quaternion `q` that rotates with angle `angle` along the y-axis of frame 1.

From the two possible solutions `q` the one is returned that is closer 
to `q_guess` (note, `q` and `-q` define the same rotation).
"""
@inline function qrot2(angle::Number; q_guess::Quaternion=NullQuaternion)::Quaternion
    q = Quaternion(0.0, sin(angle / 2), 0.0, cos(angle / 2))

    return dot(q, q_guess) >= 0 ? q : -q
end


"""
    q = ModiaMath.qrot3(angle; q_guess = NullQuaternion)

Return Quaternion `q` that rotates with angle `angle` along the z-axis of frame 1.

From the two possible solutions `q` the one is returned that is closer 
to `q_guess` (note, `q` and `-q` define the same rotation).
"""
@inline function qrot3(angle::Number; q_guess::Quaternion=NullQuaternion)::Quaternion
    q = Quaternion(0.0, 0.0, sin(angle / 2), cos(angle / 2))

    return dot(q, q_guess) >= 0 ? q : -q
end


absoluteRotation(q1::Quaternion, q_rel::Quaternion)::Quaternion = 
     (@SMatrix [ q_rel[4]  q_rel[3] -q_rel[2] q_rel[1]; 
                -q_rel[3]  q_rel[4]  q_rel[1] q_rel[2]; 
                 q_rel[2] -q_rel[1]  q_rel[4] q_rel[3]; 
                -q_rel[1] -q_rel[2] -q_rel[3] q_rel[4]]) * q1


"""
    q = ModiaMath.qrot123(angle1, angle2, angle3)

Return Quaternion `q` by rotating with angle1 along the x-axis of frame 1,
then with angle2 along the y-axis of this frame and then with angle3 along
the z-axis of this frame.

From the two possible solutions `q` the one is returned that is closer 
to `q_guess` (note, `q` and `-q` define the same rotation).
"""
qrot123(angle1::Number, angle2::Number, angle3::Number)::Quaternion = absoluteRotation(absoluteRotation(qrot1(angle1), qrot2(angle2)), qrot3(angle3))



"""
    q = ModiaMath.qrot_e(e, angle; q_guess = NullQuaternion)

Return Quaternion `q` that rotates with angle `angle` along unit axis `e`.
This function assumes that `norm(e) == 1`.

From the two possible solutions `q` the one is returned that is closer 
to `q_guess` (note, `q` and `-q` define the same rotation).
"""
@inline function qrot_e(e::Vector3D, angle::Number; q_guess::Quaternion=NullQuaternion)::Quaternion
    sa = sin(angle / 2)
    q = Quaternion(e[1] * sa, e[2] * sa, e[3] * sa, cos(angle / 2))

    return dot(q, q_guess) >= 0 ? q : -q
end
qrot_e(e::AbstractVector, angle::Number)::Quaternion = qrot_e(Vector3D(e), convert(Float64, angle))


"""
    q = ModiaMath.qrot_nxy(nx, ny)

It is assumed that the two input vectors `nx` and `ny` are resolved in frame 1 and
are directed along the x and y axis of frame 2.
The function returns the Quaternion `q` to rotate from frame 1 to frame 2. 

The function is robust in the sense that it returns always a Quaternion `q`,
even if `ny` is not orthogonal to `nx` or if one or both vectors have zero length.
This is performed in the following way: 
If `nx` and `ny` are not orthogonal to each other, first a unit vector `ey` is 
determined that is orthogonal to `nx` and is lying in the plane spanned by 
`nx` and `ny`. If `nx` and `ny` are parallel or nearly parallel to each other 
or `ny` is a vector with zero or nearly zero length, a vector `ey` is selected
arbitrarily such that `ex` and `ey` are orthogonal to each other. 
If both `nx` and `ny` are vectors with zero or nearly zero length, an
arbitrary Quaternion `q` is returned.

# Example

```julia
using Unitful
import ModiaMath

q1 = ModiaMath.qrot1(90u"°")
q2 = ModiaMath.qrot_nxy([1  , 0, 0], [0  , 0, 1  ])
q3 = ModiaMath.qrot_nxy([0.9, 0, 0], [1.1, 0, 1.1])
isapprox(q1,q2)   # returns true
isapprox(q1,q3)   # returns true
```
"""
qrot_nxy(nx, ny)::Quaternion = from_R(rot_nxy(nx, ny))




resolve1(q::Quaternion, v2::Vector3D)::Vector3D = 
    2 * ((q[4] * q[4] - 0.5) * v2 + dot(q[1:3], v2) * q[1:3] + q[4] * cross(q[1:3], v2))
           
resolve2(q::Quaternion, v1::Vector3D)::Vector3D = 
    2 * ((q[4] * q[4] - 0.5) * v1 + dot(q[1:3], v1) * q[1:3] - q[4] * cross(q[1:3], v1))

relativeRotation(q1::Quaternion, q2::Quaternion)::Quaternion = 
     (@SMatrix [ q1[4]  q1[3] -q1[2] -q1[1];
                -q1[3]  q1[4]  q1[1] -q1[2];
                 q1[2] -q1[1]  q1[4] -q1[3];
                 q1[1]  q1[2]  q1[3]  q1[4]]) * q2

inverseRotation(q::Quaternion)::Quaternion = Quaternion(-q[1], -q[2], -q[3], q[4])

