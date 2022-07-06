const oldCoords = zeros(Float64, 2)
const oldCoords2 = zeros(Float64, 2)

"""
Compute the velocity vector of a particle test at X, V of the instantaneous field
# Arguments
+ X the test particle X position
+ V the test particle V position
# Return
+ velocity vector of the test particle
"""
function instantaneous_velocity_XV(X::Float64, V::Float64, t::Float64)

    #the velocity init at zero
    velocity = zeros(2)

    #faster sin & cos
    s, c = sincos(OMEGA*t)

    #backward contribution
    velocity[1] = -V*s^2 - X*OMEGA*s*c
    velocity[2] = OMEGA^2*X*c^2 + V*OMEGA*s*c
    
    #coupled contribution
    for n = 1:NBATHPART

        velocity[1] += s * ((X - tabX[n])*c + (V-tabV[n])*s)/abs((X - tabX[n])*c + (V-tabV[n])*s) * tabm[n] * G
        velocity[2] -= c * ((X - tabX[n])*c + (V-tabV[n])*s)/abs((X - tabX[n])*c + (V-tabV[n])*s) * tabm[n] * G

    end

    return velocity

end

"""
Compute the velocity vector of a particle test at X, v of the mean field
# Arguments
+ X the test particle X position
+ V the test particle V position
# Return
+ velocity vector of the test particle
"""
function mean_velocity_XV(X::Float64, V::Float64)

    velocity = zeros(2)

    #backward contribution
    velocity[1] = -1/2*V
    velocity[2] = 1/2*OMEGA^2*X

    #coupled contribution
    for n = 1:NBATHPART
        velocity[1] += ((2*G*tabm[n])/pi)*((V - tabV[n])/(sqrt((tabX[n] - X)^2 + (1/(OMEGA^2))*(tabV[n] - V)^2)*OMEGA))
        velocity[2] -= ((2*G*tabm[n])/pi)*((X - tabX[n])/(sqrt((tabX[n] - X)^2 + (1/(OMEGA^2))*(tabV[n] - V)^2)))
    end

    return velocity

end

"""
Make the numerical derivative for computing velocity
# Arguments
+ X the test particle X position
+ V the test particle V position
# Return
+ velocity vector of the test particle
"""
function numerical_velocity(X::Float64, V::Float64)

    #classical central difference
    velocity = zeros(2)

    velocity[1] = (X - oldCoords2[1])/(2*DT)
    velocity[2] = (V - oldCoords2[2])/(2*DT)

    oldCoords2[1] = oldCoords[1]
    oldCoords2[2] = oldCoords[2]

    oldCoords[1] = X
    oldCoords[2] = V

    return velocity
end

"""
Compute the norm of the velocity
+ X the test particle X position
+ V the test particle V position
# Return
+ norm of the test particle's velocity
"""
function norm_velocity_XV(X::Float64, V::Float64, t::Float64)

    velocity = velocity_XV(X, V, t)
    return sqrt(velocity[1]^2 + velocity[2]^2)

end

"""
Compute the dot product of radial direction and velocity
+ X the test particle X position
+ V the test particle V position
# Return
+ radial component of the test particle's velocity
"""
function radial_velocity_XV(X::Float64, V::Float64)
    
    #get the radial vector
    er = zeros(2)
    er[1] = X
    er[2] = V

    #normalize it
    norm_er = sqrt.(er[1]^2+er[2]^2)

    er[1] = er[1]/norm_er
    er[2] = er[2]/norm_er

    #get the velocity
    velocity = mean_velocity_XV(X, V)

    #dot product
    radial_velocity = velocity[1]*er[1] + velocity[2]*er[2]

    return radial_velocity

end