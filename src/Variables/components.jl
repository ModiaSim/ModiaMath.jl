# License for this file: MIT (expat)
# Copyright 2017-2018, DLR Institute of System Dynamics and Control
#
# This file is part of module 
#   ModiaMath.Variables (ModiaMath/Variables/_module.jl)
#

# Definitions used to construct hierarchical variable name spaces

mutable struct ComponentInternal <: ModiaMath.AbstractComponentInternal
    name::Symbol                                                     # Name of component
    within::Union{ModiaMath.AbstractComponentWithVariables,Nothing}  # Component in which component is present (if within==nothing, object is not within a component)

    ComponentInternal(name=NoNameDefined, within=nothing) = new(Symbol(name), within)
end

isInComponent(internal::ModiaMath.AbstractComponentInternal) = typeof(internal.within) != Nothing
isNotInComponent(internal::ModiaMath.AbstractComponentInternal) = typeof(internal.within) == Nothing

isInComponent(component::ModiaMath.AbstractComponentWithVariables) = typeof(component._internal.within) != Nothing
isNotInComponent(component::ModiaMath.AbstractComponentWithVariables) = typeof(component._internal.within) == Nothing


"""
    name = ModiaMath.componentName(
               component::ModiaMath.AbstractComponentWithVariables)

Return **name** of component (without the leading path) as Symbol. 
"""
componentName(component::ModiaMath.AbstractComponentWithVariables) = component._internal.name


"""
    name = ModiaMath.fullName(component::ModiaMath.AbstractComponentWithVariables)

Return **full path name** of component (including root name) as Symbol.
"""
function fullName(component::ModiaMath.AbstractComponentWithVariables)::Symbol
    name     = component._internal.name
    internal = component._internal
    while isInComponent(internal)
        internal = internal.within._internal
        if !(internal.name ≡ NoNameDefined)
            name = Symbol(internal.name, ".", name)
        end
    end
    return name
end


"""
    ModiaMath.instanceName(component::ModiaMath.AbstractComponentWithVariables)

Return **instance name** of component (= full name but without root name) as Symbol.
"""
function instanceName(component::ModiaMath.AbstractComponentWithVariables)::Symbol
    name     = component._internal.name
    internal = component._internal
    while isInComponent(internal)
        internal = internal.within._internal
        if !(internal.name ≡ NoNameDefined)
            if isInComponent(internal)
                name = Symbol(internal.name, ".", name)
            else
                # Do not include the top-most name (because class-name, not instance-name)
                break
            end
        end
    end
    return name
end


function Base.show(io::IO, component::ModiaMath.AbstractComponentWithVariables)
    for c in fieldnames(typeof(component))
        field = getfield(component, c)
        if typeof(field) <: ModiaMath.AbstractComponentWithVariables
            println(io, "\n   ", c, " = ", field)
        elseif typeof(field) <: AbstractVector
            for i in 1:length(field)
                println(io, "   ", c, "[", i, "] = ", field[i])
            end
        elseif !( typeof(field) <: ModiaMath.AbstractComponentInternal )
            # Print fields, but not "_internal::ModiaMath.AbstractComponentInternal"
            println(io, "   ", c, " = ", field)
        end
    end

    print(io, "   )")
end



"""
    initComponent!(within::ModiaMath.AbstractComponentWithVariables, component, name)

Initialize `component` with `name` (of String or Symbol type) for component
`within`. If `component::ModiaMath.AbstractComponentWithVariables` the `within` object
and `name` are stored in `component._internal`.

Additionally, and for all other types of `component` the following statement is executed:

```julia
   setfield!(within, Symbol(name), component)
```

"""
function initComponent!(within::ModiaMath.AbstractComponentWithVariables, component, name)
    # within.name = component
    setfield!(within, Symbol(name), component)
end

function initComponent!(within::ModiaMath.AbstractComponentWithVariables, component::ModiaMath.AbstractComponentWithVariables, name)
    component._internal.within = within
    component._internal.name   = Symbol(name)

    # within.name = component
    setfield!(within, Symbol(name), component)
end



#------------------------- Macro to generate component declaration
using Base.Meta:quot, isexpr

"""
    @component ComponentName(arguments) begin ... end 

Macro to generate a `mutable struct ComponentName`. An instance of this struct
can be used as `simulationModel` in constructor [`ModiaMath.SimulationModel`](@ref)
and can then be simulated with [`ModiaMath.simulate!`](@ref).

The `arguments` must be **keyword arguments** and are used as keyword arguments for the
generated constructor function `ComponentName. The code in `begin ... end` is 
basically the body of the generated constructor function.
All left-hand-side (scalar or vector) symbols present between `begin ... end`,
as well as all keyword-arguments `arguments` are declared as fields in struct `ComponentName`.


# Examples
```julia
using ModiaMath

@component Pendulum(;L=1.0, m=1.0, d=0.1, g=9.81, phi0_deg=90.0) begin
   @assert(L > 0.0)
   @assert(m > 0.0)
   @assert(d >= 0.0)
  
   phi = RealScalar(..)
   w   = RealScalar(..)
   a   = RealScalar(..) 
   r   = RealSVector{2}(..)
end

pendulum = Pendulum(L=1.8)
```

The symbols  `L, m, d, g, phi0_deg, phi, w, a, r` are used as fields in `pendulum`
(so `pendulum.L = 1.8`).
"""
macro component(head, top_ex)
    # @component Name(arguments) begin 
    #    component1 = Component1(..)
    #    ...
    # end
    #    is mapped to:
    # mutable struct Name <: ModiaMath.AbstractComponentWithVariables
    #    <Declarations of arguments>
    #    ...
    #    Component1
    #    ...
    #    function Name(<arguments>)
    #       ##this  = new(ModiaMath.ComponentInternal(:Name, nothing), <arguments>)
    #       < equations top_ex >   
    #       ...
    #       ModiaMath.initComponent(##this, component1, :component1)
    #       ...
    #       return ##this
    #    end
    # end

    if typeof(head) == Symbol
        component_symbol = head
        fc = Expr(:call, head)
    else
        @assert isexpr(head, :call)
        component_symbol = head.args[1]::Symbol
        fc = head
    end


    # Build:  internalComponentName = ModiaMath.Component(component_name)
    # internalComponentName = gensym()
    # push!(eq, :( a = ModiaMath.Component($component_name)) )

    # Inspect top_ex find
    #   symbol = < statements >
    # and store "symbol" in "names"
    this = gensym("this")
    eqInternal = []
    eq      = []
    push!(eqInternal, :( $this = new(ModiaMath.ComponentInternal($(quot(component_symbol)))) ))
    namesSet = DataStructures.OrderedSet{Any}() # Set of left-hand side symbols
    helpArray = []

    # Add keyword arguments to namesSet
    if length(fc.args) <= 1
        # No keyword arguments
   elseif typeof(fc.args[2]) == Symbol
        # Positional arguments present, but no keyword arguments
   elseif fc.args[2].head == :parameters
        # Keyword arguments present
        kwargs = fc.args[2].args
        # println("... keyword arguments = ", kwargs)
        for kw in kwargs
            # println("... kw = ", kw, ", kw.head = ", kw.head, ", kw.args = ", kw.args, ", name = ", kw.args[1])
            # println("... kw.args[1] = ", kw.args[1], ", kw.args[2] = ", kw.args[2], ", typeof(..) = ", typeof(kw.args[2])) 
            symbol     = kw.args[1]
            # symbolType = typeof(kw.args[2])
            if in(symbol, namesSet)
                error("\nfrom @component ", component_symbol, ": keyword argument",
                  "\n   ", kw,
                  "\nis defined twice.")
            end
         
            #push!(namesSet, :( $symbol::$symbolType ) )   
            push!(namesSet, symbol)   
            push!(helpArray, :($this.$symbol = $symbol))      
        end
        # push!(fcp.args[2].args, Expr(:kw, :sceneOptions, :nothing))
    else
        error("Unknown error in components.jl: typeof(fc.args[2]) = ", typeof(fc.args[2]))
    end


    for ex in top_ex.args
        # dump(ex)
        push!(eq, ex)
        if isexpr(ex, :(=), 2) && typeof(ex.args[1]) == Symbol
            symbol = ex.args[1]
            if in(symbol, namesSet)
                error("\nfrom @component ", component_symbol, ": Left hand side symbol in",
                  "\n   ", ex,
                  "\nis defined twice.")
            end
            
            #rhs = ex.args[2]
            #if rhs.head == :call
            #   fcname = Symbol(rhs.args[1].args)
            # ??? how to figure out that fcname is a constructor of a stuct ???? If this could be figured out, fcname could be used as type in the struct (symbol::typeofSymbol)
            #end       
            push!(namesSet, symbol)
            push!(eq, :( ModiaMath.initComponent!($this, $symbol, $(quot(symbol))) ))
        end
    end
    names = collect(namesSet)

    # Build:
    #   mutable struct component_symbol <: ModiaMath.AbstractComponentWithVariables
    #      _internal::ModiaMath.ComponentInternal
    #      names[1]
    #      names[2]
    #       ...
    #      function component_symbol(...)
    #         $(eq...)
    #         ModiaMath.initComponentInternal!( new(...) )
    #      end
    #   end
    code = :(mutable struct $component_symbol <: ModiaMath.AbstractComponentWithVariables; end)
    # println("code.args = ", code.args)
    code.args[3].args = vcat(:( _internal::ModiaMath.ComponentInternal ), names)

    code2 = code.args[3].args
    push!(code2, :(function fc(); end))
    i = length(code2)
    code2[i].args[1] = fc

    # Generate constructor function
    code2[i].args[2] = quote
        $(eqInternal...)
        $(helpArray...)
        $(eq...)
        return $this
    end

    # Print generated function (without "#....\n")
    # println( replace(sprint(showcompact,code) , r"# .*\n", "\n") )
    # println( sprint(showcompact,bcode) )

    return esc(code)
end