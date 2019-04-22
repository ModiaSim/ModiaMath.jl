

function replace_nan(x)
    for i = eachindex(x)
        if isnan(x[i])
            x[i] = zero(x[i])
        end
    end
end


function time_eachindex(x)
    @time replace_nan(x)
end


x = rand(50,50,50)
[x[i] = NaN for i = (1,5,7,15,25,70)]
x1 = copy(x)
x2 = copy(x)

# compile...

time_eachindex(x)

# time

time_eachindex(x2)
