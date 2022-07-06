
"""
Compute the Gamma mean at given X and V
# Return
+ Gamma mean 
"""
function get_Gamma_mean_at(X::Float64, V::Float64)

    E_grav = 0.0
    for n=1:NBATHPART
        Xn = tabX[n]
        Vn = tabV[n]
        #Xn, Vn = xv_to_XV!(tabx[n], tabv[n], CURRENT_TIME[1])
        E_grav += (G*2.0*tabm[n]/pi)*sqrt((X - Xn)^2 + (1/OMEGA)^2*(V - Vn)^2)
    end

    E_kin = -0.25*X^2*OMEGA^2 - 0.25*V^2

    return E_grav + E_kin

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