module Runtests


@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    # Desired:
    #   using Test
    # 
    # In order that Test need not to be defined in the user environment, it is included via ModiaMath:
    using ModiaMath.Test
end

@testset "Test ModiaMath" begin
    include(joinpath("result", "_includes.jl"))
    include(joinpath("variables", "_includes.jl"))
    include(joinpath("frames", "_includes.jl"))
    include(joinpath("simulation", "_includes.jl"))
end

end