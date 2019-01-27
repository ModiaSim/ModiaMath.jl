module test_Quaternion

import ModiaMath

#  Desired:
#    using Test
#    using StaticArrays
#    using Unitful
#  
#  In order that these packages need not to be defined in the user environment, they are included via ModiaMath:
using ModiaMath.Test
using ModiaMath.StaticArrays
using ModiaMath.Unitful


angle1 = pi / 2
angle2 = 90u"°"
    
q1a = ModiaMath.qrot1(angle1)
q1b = ModiaMath.qrot1(angle2)
q1c = ModiaMath.qrot_e([1,0,0], angle1)
q1d = ModiaMath.qrot_e([1,0,0], angle2)
q1e = ModiaMath.qrot_nxy([1,0,0], [0,0,1])

q2a = ModiaMath.qrot2(angle1)
q2b = ModiaMath.qrot2(angle2)
q2c = ModiaMath.qrot_e([0,1,0], angle1)
q2d = ModiaMath.qrot_e([0,1,0], angle2)

q3a = ModiaMath.qrot3(angle1)
q3b = ModiaMath.qrot3(angle2)
q3c = ModiaMath.qrot_e([0,0,1], angle1)
q3d = ModiaMath.qrot_e([0,0,1], angle2)

q4a = ModiaMath.qrot123(angle1, angle1, angle1)
R4a = ModiaMath.rot123(angle1, angle1, angle1)

q5  = ModiaMath.qrot_nxy([1,0,0], [0,0,1])

v1  = SVector{3,Float64}(2.0, -1.0, 4.0)

v2a = ModiaMath.resolve2(q1a, v1)
v1a = ModiaMath.resolve1(q1a, v2a)

q6a = ModiaMath.absoluteRotation(q1a, q2a)
q6b = ModiaMath.qrot123(angle1, angle1, 0.0)
q7  = ModiaMath.relativeRotation(q1a, q6b)


@testset "ModiaMath.Frames: test Quaternions" begin 
    @test isapprox(q1b, q1a)
    @test isapprox(q1c, q1a)
    @test isapprox(q1d, q1a)
    @test isapprox(q1e, q1a)

    @test isapprox(q2b, q2a)
    @test isapprox(q2c, q2a)
    @test isapprox(q2d, q2a)

    @test isapprox(q3b, q3a)
    @test isapprox(q3c, q3a)
    @test isapprox(q3d, q3a)

    @test isapprox(R4a, ModiaMath.from_q(q4a))
 
    @test isapprox(q5, q1a)

    @test isapprox(v1, v1a)

    @test isapprox(q6a, q6b)
    @test isapprox(q7, q2a)

    q4d = ModiaMath.inverseRotation(q4a)
    q4e = ModiaMath.inverseRotation(q4d)
    @test isapprox(q4a, q4e)

    q4d = ModiaMath.from_R(R4a)
    R4d = ModiaMath.from_q(q4d)
    @test isapprox(R4a, R4d)
    @test isapprox(q4a, q4d)
end

end