# License for this file: MIT (expat)
# Copyright 2019, DLR Institute of System Dynamics and Control

"""
    module IdealClutch

Model of a mass controlled by a two-point controller leading to state events.
  
# Physical model

```
# Inertias
J1*der(w1) = emf.tau - clutch.tau
J2*der(w2) =           clutch.tau

# Clutch
engaged = time < 100 or time >= 300
if engaged
    0 = w2 - w1
else
    0 = clutch.tau
end

# Electric circuit
R.v   := V0 - C.v
R.i   := R.v/R
emf.v := emf.k*w1
C.v   := emf.v
C.i   := C*der(C.v)

tau1 = emf.k*emf.i
R.i  = C.i + emf.i

---------------------------
Sorted and solved equations:

input : x     = [w1     , w2     , integral(clutch_tau)]
        der_x = [der(w1), der(w2), clutch.tau]
output: residue

emf.v    := emf.k*w1
C.v      := emf.v
der(C.v) := emf.k*der(w1)
C.i      := C*der(C.v)
R.v      := V0 - C.v
R.i      := R.v/R
emf.i    := R.i - C.i
emf.tau  := emf.k*emf.i

residue1 := J1*der(w1) - (emf.tau - clutch.tau)
residue2 := J2*der(w2) - clutch.tau
residue3 := if engaged then w2-w1 else clutch.tau

is_constraint[1] := false
is_constraint[2] := false
is_constraint[3] := engaged
------

            
```

"""
module IdealClutch

import ModiaMath

const T1 = 100.0
const T2 = 300.0
		
mutable struct Model <: ModiaMath.AbstractSimulationModel
    simulationState::ModiaMath.SimulationState
   
    # Parameters   
    V0::Float64
    J1::Float64
    J2::Float64
    R::Float64
    C::Float64
    emf_k::Float64
	T_next::Float64
	engaged::Bool
      
    function Model(;V0=10.0, J1 = 0.1, J2 = 0.4, w1_start = 0.0, w2_start = 10.0, R=10.0, C=2.0, emf_k = 0.25)
        simulationState = ModiaMath.SimulationState("IdealClutch", getModelResidues!, [w1_start, w2_start, 0.0], getVariableName;
		                                            structureOfDAE = ModiaMath.DAE_LinearDerivativesAndConstraints,
													is_constraint = [false, false, true],
                                                    nz=1, nw=3)									
													
        new(simulationState, V0, J1, J2, R, C, emf_k, T1, true)
    end
end 

getVariableName(model, vcat, vindex) = ModiaMath.getVariableName(model, vcat, vindex;
                                                                 xNames    = ["inertia1.w", "inertia2.w", "integral(clutch.tau)"],
																 derxNames = ["der(inertia1.w)", "der(inertia2.w)", "clutch.tau"],
                                                                 wNames    = ["capacitor.v", "capacitor.i", "engaged"])
   
function getModelResidues!(m::Model, t::Float64, _x::Vector{Float64}, _derx::Vector{Float64}, _r::Vector{Float64}, _w::Vector{Float64})
    sim = m.simulationState
	
	if ModiaMath.isInitial(sim)
	    m.T_next  = T1
    	m.engaged = true
        ModiaMath.setNextEvent!(sim, m.T_next)
    end
		
	
    if ModiaMath.isEvent(sim)
	    # engaged = time < 100 or time >= 300

		if t >= m.T_next
            if t < T2		
		        m.T_next  = T2
			    m.engaged = false
            else
                m.T_next  = Inf
                m.engaged = true 				
	        end
	        ModiaMath.setNextEvent!(sim, m.T_next)
		end
		
		is_constraint    = ModiaMath.get_is_constraint(sim)
        is_constraint[3] = m.engaged 	
    end 	
	engaged = m.engaged


	V0    = m.V0
	J1    = m.J1
	J2    = m.J2
	R     = m.R
	C     = m.C 
	emf_k = m.emf_k
  
    w1     = _x[1]
    w2     = _x[2]
    der_w1 = _derx[1]
    der_w2 = _derx[2]
    tau    = _derx[3]

    C_v     = emf_k*w1
    der_C_v = emf_k*der_w1
    C_i     = C*der_C_v
    R_v     = V0 - C_v
    R_i     = R_v/R
    emf_i   = R_i - C_i
    emf_tau = emf_k*emf_i

    _r[1] = J1*der_w1 - (emf_tau - tau)
    _r[2] = J2*der_w2 - tau
    _r[3] = engaged ? w2-w1 : tau

    _w[1] = C_v 
	_w[2] = C_i 
	_w[3] = engaged ? 1.0 : 0.0

    #println("... time = $t, engaged = $engaged, w1 = $w1, w2 = $w2, r3 = $(_r[3])")
    return nothing
end
 
end