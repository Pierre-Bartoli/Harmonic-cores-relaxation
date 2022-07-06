"""
Perform the backward mapping from (J, theta) to (X, V)
# Arguments
+ J the action
+ theta the angle
# Return
+ (X, V) couple
"""
function JT_to_XV!(J::Float64, theta::Float64)

    #compute the sin cos
    s, c = sincos(theta)

    #compute X & V
    X = sqrt(2.0*J/OMEGA)*c
    V = sqrt(2.0*J*OMEGA)*s

    return [X, V]
end

"""
Perform the fordward mapping from (X, V) to (J, theta)
# Argument
+ X position
+ V velocities
# Return
+ (J, theta) couple with theta from 0 to 2pi
"""
function XV_to_JT!(X::Float64, V::Float64)

    #Compute action angle
    J = (X^2*OMEGA^2 + V^2)/(2.0*OMEGA)
    theta = atan(V, X) + pi

    return [J, theta]

end
