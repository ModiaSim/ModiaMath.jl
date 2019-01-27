module runtestsWithModiaAndModia3D

import Modia
import Modia3D
import ModiaMath

# Desired:
#   using Test
# 
# In order that Test need not to be defined in the user environment, it is included via ModiaMath:
using ModiaMath.Test


@testset "Test Modia and Modia3D" begin
    include("$(Modia.ModiaDir)/test/runtests.jl")
    include("$(Modia3D.path)/test/runtests.jl")

    println("\n... close all open figures.")
    ModiaMath.closeAllFigures()
end


end