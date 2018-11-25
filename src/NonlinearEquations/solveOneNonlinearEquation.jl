# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module
#   ModiaMath.NonlinearEquations (ModiaMath/NonlinearEquations/_module.jl)
#

"""
    y = solveOneNonlinearEquation(f, u_min, u_max; u_nominal=1.0, tolerance=100.0*eps())

This function determines the solution of one non-linear algebraic equation `y=f(u)`
in one unknown `u` in a reliable and efficient way. It is one of the best
numerical algorithms for this purpose. As input, the nonlinear function `f(u)`
has to be given, as well as an interval `u_min`, `u_max` that contains the solution,
that is `f(u_min)` and `f(u_max)` must have a different sign. 
The function computes iteratively smaller intervals
in which a sign change is present and terminates when the following test
is fulfilled (`u_nominal` is the nominal value of `u`, so the order of magnitude of `u`
around the solution and `tolerance` is the relative tolerance;
for example, tolerance=1e-10 means that the solution is correct up to 10 digits):

```
absTol = 0.1*u_nominal*tolerance
if abs(length_of_u_interval) <= max(tolerance*abs(u), absTol) || f(u) == 0.0
    # root u found (interval is small enough)
```

The interval reduction is performed
using inverse quadratic interpolation (interpolating with a quadratic polynomial
through the last three points and computing the zero). If this fails, bisection
is used, which always reduces the interval by a factor of two. The inverse quadratic
interpolation method has superlinear convergence. This is roughly the same convergence
rate as a globally convergent Newton method, but without the need to compute derivatives
of the non-linear function. The solver function is a direct mapping of the Algol 60 procedure
`zero` to Julia, from:

- Brent R.P. (1973): *Algorithms for Minimization without derivatives.*
  Prentice Hall, 1973, pp. 58-59.\\
  Download: [http://wwwmaths.anu.edu.au/~brent/pd/rpb011i.pdf](http://wwwmaths.anu.edu.au/~brent/pd/rpb011i.pdf)\\
  Errata and new print: [http://wwwmaths.anu.edu.au/~brent/pub/pub011.html](http://wwwmaths.anu.edu.au/~brent/pub/pub011.html)


# Examples
```julia
import ModiaMath: solveOneNonlinearEquation

fun(u; w=1.0) = 3*u - sin(w*u) - 1
u_zero = solveOneNonlinearEquation(u->fun(u; w=3.0), 0.0, 5.0)
println("residue = ", fun(u_zero; w=3.0))
```

# Remarks

- The interface was made such that it is identical to function
  [Modelica.Math.Nonlinear.solveOneNonlinearEquation](https://doc.modelica.org/Modelica%203.2.3/Resources/helpDymola/Modelica_Math_Nonlinear.html#Modelica.Math.Nonlinear.solveOneNonlinearEquation)
  in order that automatic translation of Modelica to Modia is simplified.
  However, the termination condition was changed:
  The original Brent algorithm uses an absolute tolerance for termination
  (abs(length_of_interval) <= 2*eps()*abs(u) + tolerance || f(u) == 0.0), whereas
  this was changed here to take a relative tolerance into account and use a similar
  definition of tolerances as used by Modia and Modelica,
  because easier to understand by a user (with relative tolerance and nominal value, instead of
  relative and absolute tolerances).

- Newer algorithms for the problem are presented in
  *Alefeld, Potra, Shi (1995):* [Algorithm 748: enclosing zeros of continuous functions](https://dl.acm.org/citation.cfm?id=210111).
  Here, an inverse cubic interpolation is used instead of an inverse quadratic interpolation
  as in case of the Brent algorithm. The numerical experiments in this article with
  15 test problems show that Brents algorithm needs about 5% more function evaluations
  for all the 15 test problems as the newly presented algorithm 4.2
  (in some cases Brents algorithm needs slightly less and in other cases
  slightly more function evaluations).

- Julia package [Roots.jl](https://github.com/JuliaMath/Roots.jl) provides various algorithms
  to solve the problem above but it seems that Brents algorithm is not yet included
  (as of Nov. 2018).
```
"""
function solveOneNonlinearEquation(f::Function, u_min, u_max; u_nominal=1.0, tolerance=100.0*eps())
  @assert(u_nominal > 0.0)
  @assert(tolerance > 0.0)
  @assert(u_max > u_min)
  u_min1::Float64     = convert(Float64, u_min)
  u_max1::Float64     = convert(Float64, u_max)
  u_nominal1::Float64 = convert(Float64, u_nominal)
  tolerance1::Float64 = convert(Float64, tolerance)
  absTol::Float64     = 0.1*u_nominal*tolerance1
  u::Float64 = 0.0
  a::Float64 = u_min1   # Current best minimum interval value
  b::Float64 = u_max1   # Current best maximum interval value
  c::Float64 = 0.0     # Intermediate point a <= c <= b"
  d::Float64 = 0.0
  e::Float64 = 0.0     # "b - a"
  m::Float64 = 0.0
  s::Float64 = 0.0
  p::Float64 = 0.0
  q::Float64 = 0.0
  r::Float64 = 0.0
  tol::Float64 = 0.0
  fa::Float64 = 0.0    # = f(a)
  fb::Float64 = 0.0    # = f(b)
  fc::Float64 = 0.0
  found::Bool = false


  # Check that f(u_min1) and f(u_max1) have different sign
  fa = f(u_min1)
  fb = f(u_max1)
  fc = fb
  if fa > 0.0 && fb > 0.0 || fa < 0.0 && fb < 0.0
    error("\nThe arguments u_min and u_max provided in the function call\n",
          "    solveOneNonlinearEquation(f,u_min,u_max)\n",
          "do not bracket the root of the single non-linear equation 0=f(u):\n",
          "  u_min = ", string(u_min1), "\n",
          "  u_max = ", string(u_max1), "\n",
          "  fa = f(u_min) = ", string(fa), "\n",
          "  fb = f(u_max) = ", string(fb), "\n",
          "fa and fb must have opposite sign which is not the case here.\n")
  end

  # Initialize variables
  c  = a
  fc = fa
  e  = b - a
  d  = e

  # Search
  while !found
    if abs(fc) < abs(fb)
      a  = b
      b  = c
      c  = a
      fa = fb
      fb = fc
      fc = fa
    end

    # Original formulation in Brents algorithm:
    # tol = 2*eps()*abs(b) + tolerance1
    tol = max(tolerance1*abs(b), absTol)
    m = (c - b)/2

    if abs(m) <= tol || fb == 0.0
      # root found (interval is small enough)
      found = true
      u = b
    else
      # Determine if a bisection is needed
      if abs(e) < tol || abs(fa) <= abs(fb)
        e = m
        d = e
      else
        s = fb/fa
        if a == c
          # linear interpolation
          p = 2*m*s
          q = 1 - s
        else
          # inverse quadratic interpolation
          q = fa/fc
          r = fb/fc
          p = s*(2*m*q*(q - r) - (b - a)*(r - 1))
          q = (q - 1)*(r - 1)*(s - 1)
        end

        if p > 0
          q = -q
        else
          p = -p
        end

        s = e
        e = d
        if 2*p < 3*m*q - abs(tol*q) && p < abs(0.5*s*q)
          # interpolation successful
          d = p/q
        else
          # use bi-section
          e = m
          d = e
        end
      end

      # Best guess value is defined as "a"
      a  = b
      fa = fb
      b  = b + (abs(d) > tol ? d : (m > 0 ? tol : -tol))
      fb = f(b)

      if fb > 0 && fc > 0 || fb < 0 && fc < 0
        # initialize variables
        c  = a
        fc = fa
        e  = b - a
        d  = e
      end
    end
  end

  return u
end