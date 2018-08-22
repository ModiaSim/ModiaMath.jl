# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module
#   ModiaMath.Frames (ModiaMath/Frames/_module.jl)
#

"""
    const ModiaMath.Vector3D = SVector{3,Float64}

Type of a vector in 3D space (e.g. position vector of the origin of a frame)
"""
const Vector3D = SVector{3,Float64}



"""
    const ModiaMath.ZeroVector3D = Vector3D(0.0, 0.0, 0.0)

Constant of a Vector3D where all elements are zero
"""
const ZeroVector3D = Vector3D(0.0, 0.0, 0.0)



"""
    M = ModiaMath.skew(e::AbstractVector)

Return the skew symmetric matrix `M::SMatrix{3,3,Float64,9}` of vector `e` (`length(e) = 3`)
"""
skew(e::AbstractVector) = @SMatrix([  0.0   -e[3]    e[2];
                                      e[3]   0.0    -e[1];
                                     -e[2]   e[1]    0.0 ])
