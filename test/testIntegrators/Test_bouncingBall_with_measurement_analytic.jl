module Test_bouncingBall_with_measurement_analytic

using Measurements

const g  = 9.81
const e  = 0.7 ± 0.1
const h0 = 1.0 ± 0.2

#=
der(h) = v  ; h(0) = h0
der(v) = -g ; v(0) = 0

when h <= 0 then
	v = -e*v
end when

Analytic solution of differential equation with
h(t_ev) = h_ev, v(t_ev) = v_ev:

	v = -g*(t-t_ev) + v_ev
	h = -g/2*(t-t_ev)^2 + v_ev*(t-t_ev) + h_ev

First phase (t_ev = 0, v_ev=0, h_ev=h0)
	v = -g*t
	h = -g/2*t^2 + h0

First event:
	hev1 = 0 -> t_ev1 = sqrt(2*h0/g)
	            v_ev1_before = -g*t_ev1
				v_ev1_after  = -e*v_ev1_before

Second phase (t_ev = t_ev1, v_ev = v_ev1_after, h_ev = 0)
    v = -g*(t-t_ev1) + v_ev1_after
    h = -g/2*(t-t_ev1)^2 + v_ev1_after*(t-t_ev1)

Second event:
	h_ev2 = 0 = -g/2*(t-t_ev1)^2 + vev1_after*(t-tev_1)
	          = (t-t_ev1)*(-g/2*(t-t_ev1) + v_ev1_after)

	-> t_ev2-t_ev1 = 2*v_ev1_after/g

	v_ev2_before = -g*(t_ev2-t_ev1) + v_ev1_after
	             = -g*2*v_ev1_after/g + v_ev1_after
	   			 = -v_ev1_after
	v_ev2_after  = -e*v_ev2_before
=#

t_ev1        = sqrt(2*h0/g)
h_ev1        = -g/2*t_ev1^2 + h0
v_ev1_before = -g*t_ev1
v_ev1_after  = -e*v_ev1_before

t_ev2        = t_ev1 +2*v_ev1_after/g
h_ev2        = -g/2*(t_ev2-t_ev1)^2 + v_ev1_after*(t_ev2-t_ev1)
v_ev2_before = -v_ev1_after
v_ev2_after  = -e*v_ev2_before

@show t_ev1
@show h_ev1
@show v_ev1_before
@show v_ev1_after
@show t_ev2
@show h_ev2
@show v_ev2_before
@show v_ev2_after


# Alternative that computes t_ev1 with Roots.jl:

println("\n using Roots:")

using Roots
f(t)  = -g/2*t^2 + h0
t_ev1 = find_zero(f,(0.1±0.0, 0.5±0.0))
h_ev1 = f(t_ev1)
@show t_ev1, h_ev1
end