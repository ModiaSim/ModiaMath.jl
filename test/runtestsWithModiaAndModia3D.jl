module runtestsWithModiaAndModia3D

import Modia
import Modia3D

include("$(Modia.ModiaDir)/test/runtests.jl")
include("$(Modia3D.path)/test/runtests.jl")

end