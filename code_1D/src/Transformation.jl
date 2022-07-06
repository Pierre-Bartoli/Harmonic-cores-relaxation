"""
Perform the backward mapping from (J, theta) to (X, V)
# Arguments
+ J the action
+ theta the angle
# Return
+ (X, V) couple
"""
function JT_to_XV(J::Float64, theta::Float64)

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
function XV_to_JT(X::Float64, V::Float64)

    #Compute action angle
    J = (X^2*OMEGA^2 + V^2)/(2.0*OMEGA)
    theta = atan(V, X) + pi

    return [J, theta]

end


"""
Backward mapping from (X, V) to (x, v)
# Argument
+ X postion
+ V velocity
# Return
+ (x, v) couple
"""
function XV_to_xv(X::Float64, V::Float64, t::Float64)

    #compute sin and cos
    s, c = sincos(OMEGA*t)

    #compute x & v
    x = X*c + (1.0/OMEGA)*V*s
    v = V*c - OMEGA*X*s

    return [x, v]

end

"""
Perform the forward mapping from (x, v) to (X, V)
# Argument
+ x position
+ v velocities
# Return
+ (X, V) couple
"""
function xv_to_XV(x::Float64, v::Float64, t::Float64)

    s, c  = sincos(OMEGA*t)

    X = x*c - (1.0/OMEGA)*v*s
    V = v*c + OMEGA*x*s

    return [X, V]

end
