module test_Variables

import ModiaMath

# Desired:
#   using Test
#
# In order that Test needs not to be defined in the user environment, it is included via ModiaMath:
using ModiaMath.Test


mutable struct Revolute <: ModiaMath.AbstractComponentWithVariables
    _internal::ModiaMath.ComponentInternal
    phi::ModiaMath.RealScalar
    w::ModiaMath.RealScalar
    a::ModiaMath.RealScalar
    tau::ModiaMath.RealScalar
    drive::Bool

    function Revolute(;phi0::Float64=0.0, w0::Float64=0.0, drive::Bool=false)
        this = new(ModiaMath.ComponentInternal(:Revolute))
        phi = ModiaMath.RealScalar(:phi, this, start=phi0, unit="rad",     fixed=true, info="Relative rotation angle",                     numericType=drive ? ModiaMath.WR : ModiaMath.XD_EXP)
        w   = ModiaMath.RealScalar("w", this, start=w0, unit="rad/s",   fixed=true, info="Relative angular velocity",     integral=phi, numericType=drive ? ModiaMath.WR : ModiaMath.XD_EXP,     analysis=ModiaMath.OnlyDynamicAnalysis)
        a   = ModiaMath.RealScalar("a", this, start=0.0, unit="rad/s^2",             info="Relative angular acceleration", integral=w, numericType=drive ? ModiaMath.WR : ModiaMath.DER_XD_EXP, analysis=ModiaMath.OnlyDynamicAnalysis) 
        tau = ModiaMath.RealScalar(:tau, this, start=0.0, unit="N*m",                 info="Driving torque",                              numericType=ModiaMath.WR,                                analysis=ModiaMath.QuasiStaticAndDynamicAnalysis)
        return this
    end
end

function Base.show(io::IO, rev::Revolute)
    print("Revolute(", 
         "\n    phi = ", rev.phi, 
         "\n    w   = ", rev.w, 
         "\n    a   = ", rev.a, 
         "\n    tau = ", rev.tau,
         "\n   )")
end


mutable struct Frame <: ModiaMath.AbstractComponentWithVariables
    _internal::ModiaMath.ComponentInternal
    r::ModiaMath.RealSVector3
    q::ModiaMath.RealSVector{4}

    derq::ModiaMath.RealSVector{4}
    v::ModiaMath.RealSVector3
    w::ModiaMath.RealSVector3

    a::ModiaMath.RealSVector3
    z::ModiaMath.RealSVector3   

    f::ModiaMath.RealSVector3
    t::ModiaMath.RealSVector3

    residue_w::ModiaMath.RealSVector3
    residue_f::ModiaMath.RealSVector3
    residue_t::ModiaMath.RealSVector3
    residue_q::ModiaMath.RealScalar

    drive::Bool

    function Frame(;r0=zeros(3), q0=[0,0,0,1], v0=zeros(3), w0=zeros(3), drive::Bool=false)
        this    = new(ModiaMath.ComponentInternal(:Frame))

        r         = ModiaMath.RealSVector3(:r, this, start=r0, unit="m",      fixed=true, info="Relative position",                         numericType=drive ? ModiaMath.WR : ModiaMath.XD_EXP)
        q         = ModiaMath.RealSVector{4}(:q, this, start=q0,              fixed=true, info="Relative quaternion",                       numericType=drive ? ModiaMath.WR : ModiaMath.XD_IMP)

        derq      = ModiaMath.RealSVector{4}(:derq, this,      unit="1/s",                info="der(q)",                        integral=q, numericType=drive ? ModiaMath.WR : ModiaMath.DER_XD_IMP, analysis=ModiaMath.OnlyDynamicAnalysis)
        v         = ModiaMath.RealSVector3(:v, this, start=v0, unit="m/s",    fixed=true, info="Relative velocity",             integral=r, numericType=drive ? ModiaMath.WR : ModiaMath.XD_IMP,     analysis=ModiaMath.OnlyDynamicAnalysis)
        w         = ModiaMath.RealSVector3(:w, this, start=w0, unit="rad/s",  fixed=true, info="Relative angular velocity",                 numericType=drive ? ModiaMath.WR : ModiaMath.XD_IMP,     analysis=ModiaMath.OnlyDynamicAnalysis)

        a         = ModiaMath.RealSVector3(:a, this,           unit="m/s^2",              info="Relative acceleration",         integral=v, numericType=drive ? ModiaMath.WR : ModiaMath.DER_XD_IMP, analysis=ModiaMath.OnlyDynamicAnalysis)
        z         = ModiaMath.RealSVector3(:z, this,           unit="rad/s^2",            info="Relative angular acceleration", integral=w, numericType=drive ? ModiaMath.WR : ModiaMath.DER_XD_IMP, analysis=ModiaMath.OnlyDynamicAnalysis)

        f         = ModiaMath.RealSVector3(:f, this,           unit="N",                  info="Driving force",                             numericType=ModiaMath.WR,                                analysis=ModiaMath.QuasiStaticAndDynamicAnalysis)
        t         = ModiaMath.RealSVector3(:t, this,           unit="N*m",                info="Driving torque",                            numericType=ModiaMath.WR,                                analysis=ModiaMath.QuasiStaticAndDynamicAnalysis)

        residue_w = ModiaMath.RealSVector3(:residue_w, this,                              info="Angular velocity residue",                  numericType=ModiaMath.FD_IMP,                            analysis=ModiaMath.OnlyDynamicAnalysis)
        residue_f = ModiaMath.RealSVector3(:residue_f, this,                              info="Momentum equation residue",                 numericType=ModiaMath.FD_IMP,                            analysis=ModiaMath.OnlyDynamicAnalysis)
        residue_t = ModiaMath.RealSVector3(:residue_t, this,                              info="Angular momentum equation residue",         numericType=ModiaMath.FD_IMP,                            analysis=ModiaMath.OnlyDynamicAnalysis)
        residue_q = ModiaMath.RealScalar(:residue_q, this,                                info="Quaternion constraint residue",             numericType=ModiaMath.FC,                                analysis=ModiaMath.OnlyDynamicAnalysis)

        this.drive = drive
        return this
    end
end

function Base.show(io::IO, frame::Frame)
    print("Revolute(", 
         "\n    r = ", frame.r, 
         "\n    q = ", frame.q,
         "\n    v = ", frame.v,
         "\n    w = ", frame.w,
         "\n    a = ", frame.a,
         "\n    z = ", frame.z,
         "\n    f = ", frame.f,
         "\n    t = ", frame.t,
         "\n   )")
end

ModiaMath.@component Robot(;phi10=0.0, phi20=0.0, var10=0.0, r0=zeros(3), q0=zeros(4)) begin
    rev1  = Revolute(phi0=phi10)
    rev2  = Revolute(phi0=phi20)
    var1  = ModiaMath.RealScalar(start=var10, numericType=ModiaMath.XA)
    res1  = ModiaMath.RealScalar(numericType=ModiaMath.FD_IMP)
    frame = Frame(r0=r0, q0=q0)
end


ModiaMath.@component Robot2(;phi10=0.0, phi20=0.0, r0=zeros(3), q0=[0.0, 0.0, 0.0, 1.0]) begin
    rev1  = Revolute(phi0=phi10, drive=true)
    rev2  = Revolute(phi0=phi20, drive=true)
    frame = Frame(r0=r0, q0=q0,  drive=true)
end


ModiaMath.@component Robot3(;phi10=0.0, phi20=0.0, phi30=-2.0) begin
    rev1 = Revolute(phi0=phi10, drive=true)
    rev2 = Revolute(phi0=phi20, drive=true)
    rev3 = Revolute(phi0=phi30)
    res1 = ModiaMath.RealScalar(numericType=ModiaMath.FC)
end



@testset "ModiaMath: test_Variables.jl" begin 

   
    @testset "Dynamic analysis" begin # ------------------------------ Dynamic analysis 
        r0 = [1.0,2.0,3.0]
        q0 = [0.5,0.5,0.0,sqrt(0.5^2 + 0.5^2)]
        robot = Robot(phi10=1.0, phi20=2.0, var10=3.0, r0=r0, q0=q0)
        robot.rev1.a.value = 10 * 2.22
        robot.rev2.a.value = 10 * 4.44
   
        println("\n... robot = ", robot) 
   
        println("\n... Print variables of robot")
        m = ModiaMath.ModelVariables(robot)
        ModiaMath.print_ModelVariables(m)
   
        println("\n... Copy start values to x")
        x         = zeros(5 + 7 + 6)
        x_fixed   = fill(false, 5 + 7 + 6)
        x_nominal = fill(2.0, 5 + 7 + 6)
        ModiaMath.copy_start_to_x!(m, x, x_fixed, x_nominal)
        x0 = [1.0, 0.0, 2.0, 0.0, 1.0, 2.0, 3.0, 0.5, 0.5, 0.0, sqrt(0.5^2 + 0.5^2), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0 ]
        x0_fixed = fill(true, 17)
        append!(x0_fixed, false)
        @test isapprox(x, x0)
        @test x_fixed == x0_fixed
  
        println("\n... Copy x and der_x to variables")
        x    = [1.11, 2.22 , 3.33, 4.44, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 1.1, 2.2, 3.3, 4.4, 5.5, 6.6 , 5.55]
        derx = [2.22, 22.2, 4.44, 44.4, 1.1, 2.2, 3.3, 4.44, 5.55, 6.66, 7.77, 1.11, 2.22, 3.33, 4.44, 5.55, 6.66, 0.0]  
        ModiaMath.copy_x_and_derx_to_variables!(0.5, x, derx, m)
        @test isapprox(x, [robot.rev1.phi.value, robot.rev1.w.value, robot.rev2.phi.value, robot.rev2.w.value,
                            robot.frame.r.value[1], robot.frame.r.value[2], robot.frame.r.value[3],
                            robot.frame.q.value[1], robot.frame.q.value[2], robot.frame.q.value[3], robot.frame.q.value[4], 
                            robot.frame.v.value[1], robot.frame.v.value[2], robot.frame.v.value[3],
                            robot.frame.w.value[1], robot.frame.w.value[2], robot.frame.w.value[3],
                            robot.var1.value], rtol=1e-15)
        @test isapprox(derx[2], robot.rev1.a.value, rtol=1e-15)
        @test isapprox(derx[4], robot.rev2.a.value, rtol=1e-15)
   
        println("\n... Copy variables to residues")
        residues = zeros(5 + 7 + 6)
        ModiaMath.copy_variables_to_residue!(m, x, derx, residues)
        println("residue = ", residues)
        @test isapprox(residues, zeros(5 + 7 + 6), atol=1e-12)
    end



    @testset "Kinematic analysis 1" begin   # ---------------------------- Kinematic analysis 1
        robot2 = Robot2(phi10=1.0, phi20=2.0)
        println("\n... robot2 = ", robot2) 
   
        println("\n... Print variables of robot2")
        m = ModiaMath.ModelVariables(robot2, analysis=ModiaMath.KinematicAnalysis)
        ModiaMath.print_ModelVariables(m)
   
        println("\n... Copy start values to x")
        x         = [1.11]
        x_fixed   = [false]
        x_nominal = [3.0]
        ModiaMath.copy_start_to_x!(m, x, x_fixed, x_nominal)
        @test isapprox(x, [0.0])
        @test x_fixed == [true]
   
        println("\n... Copy x and der_x to variables")
        x    = fill(2.1, 1)
        derx = fill(3.2, 1)  
        ModiaMath.copy_x_and_derx_to_variables!(0.5, x, derx, m)
        @test isapprox(x, [m.var[2].value], rtol=1e-15)    # var[2] = _dummy_x
   
        println("\n... Copy variables to residues")
        residues = zeros(m.nx)
        ModiaMath.copy_variables_to_residue!(m, x, derx, residues)
        @test isapprox(residues, [3.2 - (-2.1)], atol=1e-12)
    end



    @testset "Kinematic analysis 2" begin   # ---------------------------- Kinematic analysis 2
        robot3 = Robot3(phi10=1.0, phi20=2.0, phi30=-2.0)
        println("\n... robot3 = ", robot3) 
   
        println("\n... Print variables of robot3")
        m = ModiaMath.ModelVariables(robot3, analysis=ModiaMath.KinematicAnalysis)
        ModiaMath.print_ModelVariables(m)
   
        println("\n... Copy start values to x")
        x = zeros(m.nx)
        x_fixed   = fill(false, m.nx)
        x_nominal = fill(10.0 , m.nx)
        ModiaMath.copy_start_to_x!(m, x, x_fixed, x_nominal)
        @test isapprox(x, [-2.0])
        @test x_fixed == [true]
   
        println("\n... Copy x and der_x to variables")
        x    = [1.11]
        derx = [2.22]
        ModiaMath.copy_x_and_derx_to_variables!(0.5, x, derx, m)
        @test isapprox(x, [robot3.rev3.phi.value], rtol=1e-15)  
   
        println("\n... Copy variables to residues")
        residues = zeros(m.nx)
        ModiaMath.copy_variables_to_residue!(m, x, derx, residues)
        @test isapprox(residues, zeros(1), atol=1e-12)
    end

end

end