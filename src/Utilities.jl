# License for module ModiaMath.Logging: MIT
# Copyright 2017-2018, DLR Institute of System Dynamics and Control


"""
    module ModiaMath.Utilities

Utility functions used in ModiaMath.

# Main developer
Martin Otter, [DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)
"""
module Utilities

export printobj


"""
    str = displayAsString(x)

Return the result of the Julia built-in command display(x) as string 
instead of displaying it on the terminal.
"""
displayAsString(x) = sprint((io, x) -> display(TextDisplay(io), x), x)


"""
    str = displayAsString2(x)

Return the result of the Julia built-in command display(x) as string 
instead of displaying it on the terminal, remove the preceding type
information and shift the whole output to the right
"""
function displayAsString2(obj)
    str  = displayAsString(obj)
    str2 = str[searchindex(str, "\n") + 1:end]
    str3 = "     " * replace(str2, "\n", "\n     ")
end


"""
    ModiaMath.printobj(name,obj)

Pretty print obj as "<name> = <display(obj)>". 
"""
function printobj(name, obj)
    str = displayAsString2(obj)
    println("\n  ", name, " =\n", str)
end

end