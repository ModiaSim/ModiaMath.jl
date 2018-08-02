module test_Interpolation

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end


using Unitful
import ModiaMath


r = [ ModiaMath.Vector3D(1,0,0),
      ModiaMath.Vector3D(0,1,0),
      ModiaMath.Vector3D(0,0,1) ]

q = [ ModiaMath.NullQuaternion,
      ModiaMath.qrot1(45u"°"),
      ModiaMath.qrot2(60u"°")  ]

path     = ModiaMath.Path(r,q)
t_end    = ModiaMath.t_pathEnd(path)
println("t_end = ", t_end)
println("path.t = ", path.t)
dt       = 0.1
stopTime = 2.0 
time     = 0.0

@testset "ModiaMath.Frame: test interpolation" begin 
   dist1 = norm(r[2] - r[1])
   dist2 = norm(r[3] - r[2])

   @test isapprox(ModiaMath.t_pathEnd(path), dist1+dist2)

   (rt2, qt2) = ModiaMath.interpolate(path, dist1)
   @test isapprox(r[2], rt2)
   @test isapprox(q[2], qt2)

   rt3 = ModiaMath.interpolate_r(path, dist1/2.0)
   @test isapprox(rt3, r[1] + (r[2] - r[1])/2 )
end


while time <= stopTime + 1e-7
   (rt, qt) = ModiaMath.interpolate(path, time*t_end/stopTime)
   println("... time = ", time, ", rt = ", rt)
   time += dt
end

end