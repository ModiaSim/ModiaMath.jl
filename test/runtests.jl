module Runtests

import ModiaMath

# Desired:
#   using Test
# 
# In order that Test need not to be defined in the user environment, it is included via ModiaMath:
using ModiaMath.Test


@testset "Test ModiaMath" begin
    include(joinpath("result", "_includes.jl"))
    include(joinpath("variables", "_includes.jl"))
    include(joinpath("frames", "_includes.jl"))
    include(joinpath("nonlinearEquations", "_includes.jl"))
    include(joinpath("simulation", "_includes.jl"))

    println("\n... close all open figures.")
    ModiaMath.closeAllFigures()
end

end