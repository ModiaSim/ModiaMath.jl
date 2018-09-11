module Runtests

import ModiaMath
using Test

@testset "Test ModiaMath" begin
    include(joinpath("result", "_includes.jl"))
    include(joinpath("variables", "_includes.jl"))
    include(joinpath("frames", "_includes.jl"))
    include(joinpath("simulation", "_includes.jl"))
end

end
