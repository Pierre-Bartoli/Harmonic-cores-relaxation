"""cumulative tabm*tabx from left to right
"""
const S1 = zeros(Float64, NPART + 1)

"""cumulative tabm from left to right
"""
const Sm1 = zeros(Float64, NPART + 1)

"""cumulative tabm*tabx from right to left
"""
const S2 = zeros(Float64, NPART + 1)

"""cumulative tabm form right to left
"""
const Sm2 = zeros(Float64, NPART + 1)

"""
Find the nearest index of x in the array by dichotomy
# Argument 
+ x element
+ tab to where to find the nearest index of x
+ tab_perm the array that gives the sorted index
# Return
+ index where x is the minimum higher position
"""
function find_index(tab_perm::Array{Int64, 1}, tab::Array{Float64, 1}, x::Float64)
    
    a = 1 #initialize min
    b = NPART #initialize max
    found = false # if the good index is fonded
    c = trunc(Int, (a + b)/2) #get the middle

    ind1 = tab_perm[1] # if x is under all particles
    if (x < tab[ind1])
        c = 1
        found = true
    end
    indN = tab_perm[NPART] # if x is upper all particles
    if (x > tab[indN])
        c = NPART + 1
        found = true
    end

    while (!found && a <= b) # dichotomy algorithm 
        indc = tab_perm[c]
        indc1 = tab_perm[c+1]
        if (tab[indc] <= x && x <= tab[indc1]) # found when x is between two indexes
            found = true
        else
            if (x > tab[indc])
                a = c + 1
            else
                b = c - 1 
            end
            c = trunc(Int, (a + b)/2)
        end
    end

    return c
end

"""
Get the Gamma at given x, and cumulative sums
# Arguments
+ S1 sum of cumulative mass * x
+ S2 reverse sum of cumulative mass * x (in reverse side)
+ Sm1 cumulative mass
+ Sm2 cumulative mass in reverse
+ x position to evalutate the Gamma
# Return
+ Gamma energy = E_kin + E_grav \n 
E_kin = -1/2x^2*omega^2 \n
E_grav = sum(abs(x - x_i))
"""
function get_Gamma_inst_at(ind::Int, x::Float64)

    #computed gravitationnal energy
    E_grav = S2[ind] - S1[ind] + x*(Sm1[ind] - Sm2[ind])

    #energy kinetic in big variables
    E_kin = -0.5*OMEGA^2*x^2
    
    return E_kin + E_grav

end

"""
Compute the Gamma mean at given X and V
# Return
+ Gamma mean 
"""
function get_Gamma_mean_at(X::Float64, V::Float64)

    #compute gravitationnal energy
    E_grav = 0.0
    for n=1:NBATHPART
        Xn = tabX[n]
        Vn = tabV[n]
        #Xn, Vn = xv_to_XV!(tabx[n], tabv[n], CURRENT_TIME[1])
        E_grav += (G*2.0*tabm[n]/pi)*sqrt((X - Xn)^2 + (1/OMEGA)^2*(V - Vn)^2)
    end

    #kinetic energy
    E_kin = -0.25*X^2*OMEGA^2 - 0.25*V^2

    return E_grav + E_kin

end

"""
Fill the cumulatives sums to compute faster the Gamma instantaneous
"""
function fill_cumulative_sums!()

    fill!(S1, 0.0)
    fill!(Sm1, 0.0)
    fill!(S2, 0.0)
    fill!(Sm2, 0.0)

    for n=1:NPART #over the NPART
        ind = tabperm[n] # permutation index at n
        indN1 = tabperm[NPART - n + 1] # because S2 is reverse take the index at N - n + 1
        S1[n + 1] = S1[n] + tabm[ind]*tabx[ind] # cumulative mass*x
        Sm1[n + 1] = Sm1[n] + tabm[ind] # Cumulative mass
        S2[NPART - n + 1] = S2[NPART - n + 2] + tabm[indN1]*tabx[indN1] # cumulative mass*x in reverse (from right to left)
        Sm2[NPART - n + 1] = Sm2[NPART - n + 2] + tabm[indN1] # cumulative mass from right to left
    end

end

"""
Compute the instantaneous gamma for given J and theta
# Arguments
+ J the action
+ theta the angle
"""
function Gamma_inst(J::Float64, theta::Float64)

    X, V = JT_to_XV(J, theta) # transformations
    x, v = XV_to_xv(X, V, get_time!()) # transformations
            
    ind = find_index(tabperm, tabx, x) # find the index
    return get_Gamma_inst_at(ind, x) #compute Gamma instantaneous

end

""" 
Compute Gamma mean for given J and theta
# Arguments
+ J the action
+ theta the angle
"""
function Gamma_mean(J::Float64, theta::Float64)

    X, V = JT_to_XV(J, theta) # transformation
    return get_Gamma_mean_at(X, V)# compute Gamma mean

end

""" 
Compute Gamma mean for given X and V
# Arguments
+ X position
+ V velocity
"""
function Gamma_mean_XV(X::Float64, V::Float64)

    return get_Gamma_mean_at(X, V)# compute Gamma mean

end