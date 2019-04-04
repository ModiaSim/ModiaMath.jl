# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control

"""
    module PendulumDAE

DAE-model of a mass point attached via a rod to a revolute joint (GGL formulation).
    
Starting equations (der(xx) = dxx/dt):

     der(x) = vx
     der(y) = vy
  m*der(vx) = -x*lambda
  m*der(vy) = -y*lambda-m*g
  x*x + y*y = L*L

Gear/Gupta/Leimkuhler equations:
   # lambda = der(lambda_int)
   # mue    = der(mue_int)
   # g      = x*x + y*y - L*L
   # G      = dg/d[x;y] = [2*x 2*y]
   
   der([x,y]) = [vx,vy] - G'*mue
   0 = m*der(vx) + lambda*x
   0 = m*der(vy) + lambda*y + m*g
   0 = x*x + y*y - L*L
   0 = 2*x*vx + 2*y*vy
   
Arguments of getResidues! function:
   _x    = [x; y; vx; vy; lambda_int; mue_int]
   _derx = [der(x); der(y); der(vx); der(vy); lambda; mue]

Modelica model: ModelicaReferenceModels.ODAEs.Pendulum   
"""
module PendulumDAE

import ModiaMath

mutable struct Model <: ModiaMath.AbstractSimulationModel
    simulationState::ModiaMath.SimulationState
   
    # Parameters
    L::Float64
    m::Float64
    g::Float64

    function Model(;L=1.0, m=1.0, g=9.81, x0=L / 2.0, y0=-0.5, x_fixed=false)   #y0=-sqrt(L*L - x0*x0))
        @assert(L > 0.0)
        @assert(m > 0.0)
        @assert(-L <= x0 <= L)
        @assert(-L <= y0 <= L)      
        simulationState = ModiaMath.SimulationState("PendulumDAE", getModelResidues!, [x0,y0,1.0,1.0,0.0,0.0], getVariableName;
                                x_fixed=[x_fixed, false, x_fixed, false, false, false],
                                is_constraint = [false,false,false,false,true,true],
								has_constraintDerivatives = true)
        new(simulationState, L, m, g)
    end
end 

getVariableName(model, vcat, vindex) = ModiaMath.getVariableName(model, vcat, vindex;
                                                         xNames=["x", "y", "vx", "vy", "lambda_int", "mue_int"],
                                                         derxNames=["der(x)", "der(y)", "der(vx)", "der(vy)", "lambda", "mue"]) 
                   
function getModelResidues!(m::Model, t::Float64, _x::Vector{Float64}, _derx::Vector{Float64}, _r::Vector{Float64}, _w::Vector{Float64})   
    x      = _x[1]
    y      = _x[2]
    vx     = _x[3]
    vy     = _x[4]
    dervx  = _derx[3]
    dervy  = _derx[4]
    lambda = _derx[5]
    mue    = _derx[6]

    _r[1] = _derx[1] - vx + x * mue
    _r[2] = _derx[2] - vy + y * mue
    _r[3] = m.m * dervx + lambda * x
    _r[4] = m.m * dervy + lambda * y + m.m * m.g
	
	if ModiaMath.compute_der_fc(m)
		_r[5] = x * _derx[1] + y * _derx[2]
		_r[6] = x * dervx + y * dervy + vx*vx + vy*vy
	else
		_r[5] = (x * x + y * y - m.L * m.L) / 2.0
		_r[6] = x * vx + y * vy
    end
    return nothing
end

end