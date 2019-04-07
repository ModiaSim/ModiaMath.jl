# License for this file: MIT (expat)
# Copyright 2019, DLR Institute of System Dynamics and Control

"""
    module Simulate_PendulumDAE

Demonstrate how a pendulum described as index 3 DAE can be defined, simulated
and plotted with ModiaMath.

The model is a mass point attached via a rod to a revolute joint where "x" and "y"
coordinates are used as describing variables. The index 3 model is transformed to
an index 1 model using the Gear-Gupta-Leimkuhler transformation:
    
Starting equations (der(v) = dv/dt):

     der(x) = vx
     der(y) = vy
  m*der(vx) = -x*lambda
  m*der(vy) = -y*lambda-m*g
  x*x + y*y = L*L

Gear/Gupta/Leimkuhler transformation:
   # lambda = der(lambda_int)
   # mue    = der(mue_int)
   # g      = x*x + y*y - L*L
   # G      = dg/d[x;y] = [2*x 2*y]
   
   der([x,y]) = [vx,vy] - G'*mue
   0 = m*der(vx) + lambda*x
   0 = m*der(vy) + lambda*y + m*g
   0 = x*x + y*y - L*L
   0 = 2*x*vx + 2*y*vy
   
Modelica model: ModelicaReferenceModels.ODAEs.Pendulum   
"""
module Simulate_PendulumDAE

using ModiaMath

@component PendulumDAE(;L=1.0, m=1.0, g=9.81, x_start=L/2.0, y_start=0.5, x_fixed=false) begin    # y_start = -sqrt(L*L - x_start*x_start)
    @assert(L > 0.0)
    @assert(m > 0.0)
  
    x      = RealScalar(           numericType=ModiaMath.XD_EXP,                  start=x_start, unit="m"    , fixed=x_fixed, nominal=1.0, info="Absolute x-position of mass")
    y      = RealScalar(           numericType=ModiaMath.XD_EXP,                  start=y_start, unit="m"    , fixed=false  , nominal=1.0, info="Absolute y-position of mass")
    der_x  = RealScalar("der(x)",  numericType=ModiaMath.DER_XD_EXP, integral=x,                 unit="m/s"  ,                             info="= der(x)")
    der_y  = RealScalar("der(y)",  numericType=ModiaMath.DER_XD_EXP, integral=y,                 unit="m/s"  ,                             info="= der(y)")
	
    vx     = RealScalar(           numericType=ModiaMath.XD_EXP,                  start=0.0,     unit="m/s"  , fixed=x_fixed, nominal=1.0, info="Absolute x-velocity of mass")
    vy     = RealScalar(           numericType=ModiaMath.XD_EXP,                  start=0.0,     unit="m/s"  , fixed=false  , nominal=1.0, info="Absolute y-velocity of mass")
    ax     = RealScalar("der(vx)", numericType=ModiaMath.DER_XD_EXP, integral=vx,                unit="m/s^2",                             info="= der(vx)")
    ay     = RealScalar("der(vy)", numericType=ModiaMath.DER_XD_EXP, integral=vy,                unit="m/s^2",                             info="= der(vy)")

    lambda        = RealScalar(numericType=ModiaMath.LAMBDA, unit="N"    , info="Cut-force in the rod")
    mue           = RealScalar(numericType=ModiaMath.MUE,    unit="1/s"  , info="Mue stabilizer (should be zero)")	
    residue_xy    = RealScalar(numericType=ModiaMath.FC,     unit="m^2"  , info="Residue of position constraint")
    residue_vx_vy = RealScalar(numericType=ModiaMath.FC,     unit="m^2/s", info="Residue of velocity constraint")
end
          

function ModiaMath.computeVariables!(p::PendulumDAE, sim::ModiaMath.SimulationState)
    x      = p.x.value
    y      = p.y.value
    vx     = p.vx.value
    vy     = p.vy.value
    lambda = p.lambda.value
    mue    = p.mue.value

    p.der_x.value = vx + x * mue
	p.der_y.value = vy + y * mue
    p.ax.value    = lambda * x/p.m	
    p.ay.value    = lambda * y/p.m + p.g

	
#	if ModiaMath.compute_der_fc(m)  
#		_r[5] = x * _derx[1] + y * _derx[2]
#		_r[6] = x * dervx + y * dervy + vx*vx + vy*vy
#	else
	p.residue_xy.value    = (x * x + y * y - p.L * p.L) / 2.0
	p.residue_vx_vy.value = x * vx + y * vy
#    end
    return nothing
end

simulationModel = ModiaMath.SimulationModel(PendulumDAE(L=0.8, m=0.5); stopTime=5.0, tolerance=1e-6)
ModiaMath.print_ModelVariables(simulationModel)
result          = ModiaMath.simulate!(simulationModel, log=true)
ModiaMath.plot(result, [("x", "y"), ("vx", "vy"), "lambda", "mue"])

end