module test_SolveOneNonlinearEquation

import ModiaMath

#  Desired:
#    using Test
#    using LinearAlgebra
#    using Unitful
#  
#  In order that these packages need not to be defined in the user environment, they are included via ModiaMath:
using ModiaMath.Test
@eval using Printf
using ModiaMath: solveOneNonlinearEquation


@testset "ModiaMath.NonlinearEquations: test solveOneNonlinearEquation" begin 
    # Test from Modelica.Math.Nonlinear.Examples.solveNonlinearEquations1
    y          = zeros(3)
    y_analytic = zeros(3)

    fun1(u)       = u^2 - 1
    y[1]          = solveOneNonlinearEquation(fun1, -0.5, 10.0)
    y_analytic[1] = 1.0

    fun2(u;w=1.0) = 3*u - sin(w*u) - 1
    y[2]          = solveOneNonlinearEquation(u->fun2(u; w=3.0), 0.0, 5.0)
    y_analytic[2] = 0.6448544035840080891877538

    fun3(u; p=zeros(2), m=0.0) = p[1] + log(p[2]*u) - m*u
    y[3]          = solveOneNonlinearEquation(u->fun3(u; p=[5.0,1.0], m=1.0), 1.0, 8.0)
    y_analytic[3] = 6.9368474072202187221643182


    for i in 1:length(y)
        println("\n... Results of Solve_SingleNonlinearEquations:")
        println("fun", string(i), ":")
        @printf("   analytical zero     = %20.16e\n", y_analytic[i])
        @printf("   numerical zero      = %20.16e\n", y[i]) 
        @printf("   absolute difference = %20.16e\n", y_analytic[i] - y[i])
        @test isapprox(y[i], y_analytic[i],  atol=1e-12)
    end


    # Test error case
    errorOccured = false
    try
       yy = solveOneNonlinearEquation(u->fun2(u; w=3.0), 0.0, 0.5)
    catch
       errorOccured = true
    end
    @test errorOccured
 

    # Test tolerance
    yy = solveOneNonlinearEquation(fun1, -0.5, 10.0; tolerance=1e-4)
    @test isapprox(yy, 1.0,  atol=1e-4)



    # Tests from https://github.com/JuliaMath/Roots.jl/test/test_find_zero.jl
    fns = [(x -> x^2 - exp(x) - 3x + 2,   0.257530285439860, -1, 0.5),
           (x -> x^3 + 4x^2 -10,          1.365230013414097, 0.5, 2.0),
           (x -> exp(x)*sin(x) - 2x - 5, -2.523245230732555, -3.0, -1.0),
           (x -> log(x) + sqrt(x) - 5,    8.309432694231571,  7, 9),
           (x -> sqrt(x) - 1/x - 3,       9.633595562832696,  8, 10),
           (x -> exp(-x) - cos(x), 1.292695719373398, 1.0, 2.0),
           (x -> cos(x)^2 - x/5, 1.085982678007472, 0.5, 2.0),
           (x -> x^10 - x^3 - x - 1, -0.674177935277052, -1.0, 0.0),
           (x -> sin(x) - x + 1, 1.934563210752024, 1.0, 20.0),
           (x -> exp(-x^2 + x + 2) - cos(x) + x^3 + 1, -1.0899423334391694, -2, 0),
           (x -> asin(x^2 -1) - x/2 + 1, 0.594810968398369, 0.3, 0.9),
           (x -> tanh(x) - tan(x), 7.068582745628732, 5.5, 7.1)
    ]

    for (i, (f, u_exact, u_min, u_max)) in enumerate(fns)
        u = solveOneNonlinearEquation(f, u_min, u_max)
        @test isapprox(u, u_exact,  atol=1e-12)
    end


    # Some more tests from https://github.com/JuliaMath/Roots.jl/test/test_find_zero.jl
    # based on http://ir.igsnrr.ac.cn/bitstream/311030/8840/1/%E4%BE%AF%E9%BA%9F%E7%A7%91(SCI)2.pdf
    # which have multiple roots
    multiplicity_tests = [
         (x -> (x - sqrt(5))^4 / ((x-1)^2 + 2),    2.236067977499790, 0, 3)
         (x -> (sin(x)^2 - 2x + 1)^5,              0.71483582544138924, 0, 1)
         (x -> (8x*exp(-x^2) -2x - 3)^8,           -1.7903531791589544, -2, -1)
         (x -> (2x*cos(x) + x^2 - 3)^10/(x^2 + 1),  2.9806452794385368, 0, 3)
         (x -> (exp(-x^2 + x + 3) - x + 2)^9,       2.4905398276083051, 0, 4)

         (x -> (exp(-x) + 2sin(x))^4,              3.1627488709263654, 2, 4)
         (x -> (log(x^2 + 3x + 5) - 2x + 7)^8,     5.4690123359101421, 3, 6)
         (x -> (sqrt(x^2 + 2x + 5) - 2sin(x) - x^2 + 3)^5,   2.3319676558839640, 0, 4)
         (x -> (x-2)^4 / ( (x-1)^2 + 1),           2.0000000000000000, 0, 4)
         (x -> abs(x - 2.5)^(15/4) * exp(x),       2.5, 0, 4)

         (x -> (sqrt(x) - 1/x - 1)^7,              2.147899035704787, 0, 4)
         (x -> (log(x) + sqrt(x) - 5)^3,           8.309432694231572, 0, 10)
         (x -> (sin(x)*cos(x) - x^3 + 1)^9,        1.117078770687451, 0, 2)
         (x -> ((x-3)*exp(x))^5,                   3.0000000000000000, 0, 5)
         (x -> (log(x) + sqrt(x^4 + 1) -2)^7,      1.222813963628973, 1, 3)
         ]

    for (i, (f, u_exact, u_min, u_max)) in enumerate(fns)
        u = solveOneNonlinearEquation(f, u_min, u_max)
        @test isapprox(u, u_exact,  atol=1e-12)
    end
end

end