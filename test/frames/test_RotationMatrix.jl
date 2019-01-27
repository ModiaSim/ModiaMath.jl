module test_RotationMatrix

import ModiaMath

#  Desired:
#    using Test
#    using LinearAlgebra
#    using StaticArrays
#    using Unitful
#  
#  In order that these packages need not to be defined in the user environment, they are included via ModiaMath:
using ModiaMath.Test
using ModiaMath.LinearAlgebra
using ModiaMath.StaticArrays
using ModiaMath.Unitful


angle1 = pi / 2
angle2 = 90u"°"
    
R1a = ModiaMath.rot1(angle1)
R1b = ModiaMath.rot1(angle2)
R1c = ModiaMath.rot_e([1,0,0], angle1)
R1d = ModiaMath.rot_e([1,0,0], angle2)
R1e = ModiaMath.rot_nxy([1,0,0], [0,0,1])

R2a = ModiaMath.rot2(angle1)
R2b = ModiaMath.rot2(angle2)
R2c = ModiaMath.rot_e([0,1,0], angle1)
R2d = ModiaMath.rot_e([0,1,0], angle2)

R3a = ModiaMath.rot3(angle1)
R3b = ModiaMath.rot3(angle2)
R3c = ModiaMath.rot_e([0,0,1], angle1)
R3d = ModiaMath.rot_e([0,0,1], angle2)

R4a = ModiaMath.rot123(angle1, angle1, angle1)
R4b = ModiaMath.rot123(angle2, angle2, angle2)
R4c = ModiaMath.rot_e([1,0,0], -angle1) *
      ModiaMath.rot_e([0,1,0], -angle1) *
      ModiaMath.rot_e([0,0,1], -angle1)

R5  = ModiaMath.rot_nxy([1,0,0], [0,0,1])

v1  = SVector{3,Float64}(2.0, -1.0, 4.0)

v2a = ModiaMath.resolve2(R1a, v1)
v1a = ModiaMath.resolve1(R1a, v2a)

R6a = ModiaMath.absoluteRotation(R1a, R2a)
R6b = ModiaMath.rot123(angle1, angle1, 0.0)
R7  = ModiaMath.relativeRotation(R1a, R6b)

@testset "ModiaMath.Frames: test RotationMatrix" begin 
    @test isapprox(R1b, R1a)
    @test isapprox(R1c, R1a)
    @test isapprox(R1d, R1a)
    @test isapprox(R1e, R1a)

    @test isapprox(R2b, R2a)
    @test isapprox(R2c, R2a)
    @test isapprox(R2d, R2a)

    @test isapprox(R3b, R3a)
    @test isapprox(R3c, R3a)
    @test isapprox(R3d, R3a)

    @test isapprox(R4b, R4a)
    @test isapprox(R4c, R4a)
 
    @test isapprox(R5, R1a)

    @test isapprox(v1, v1a)

    @test isapprox(R6a, R6b)
    @test isapprox(R7, R2a)

    R4d = ModiaMath.inverseRotation(R4a)
    R4e = ModiaMath.inverseRotation(R4d)
    @test isapprox(R4a, R4e)

    angle4 = 45u"°"
    e      = normalize([1.0, 1.0, 1.0])
    R      = ModiaMath.rot_e(e, angle4)
    v3     = [1.0, 2.0, 3.0]
    v4     = ModiaMath.resolve2(R, v3)
    angle5 = ModiaMath.planarRotationAngle(e, v3, v4)
    isapprox(angle4, angle5)

    isapprox(-ModiaMath.eAxis(3), ModiaMath.eAxis(-3))
end

end