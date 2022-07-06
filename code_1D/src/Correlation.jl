"""
This function compute correlation
# Argument 
+ Gamma_mean_1 the first value of the correlation
+ Gamma_mean_2 the second value of the correlation
"""
function correlation(Gamma_mean_1::Float64, Gamma_mean_2::Float64)

   return (Gamma_mean_1)*(Gamma_mean_2)

end

"""
This function compute the overall average of the correlation
# return
correlation
"""
function correlation_mean(correlation::Array{Float64, 1})

    #make the average
    correlation_mean = 0 
    for iseed=1:NSEED
        correlation_mean += correlation[iseed]
    end

    correlation_mean /= NSEED
    return correlation_mean

end