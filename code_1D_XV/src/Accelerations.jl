"""
Compute the velocities of X and P of the particles
"""
function tab_dot_mean!()

    for i=1:NPART
        tabXdot[i] = -0.5*OMEGA*tabV[i]
        tabVdot[i] =  0.5*OMEGA^2*tabX[i]
    end
    for i = 1:NPART
        for j = (i+1):NPART
            inv_norm = (2.0*G/pi)/sqrt((tabX[i] - tabX[j])^2 + ((1.0/OMEGA)*(tabV[i] - tabV[j]))^2)
            Xdot_ij = 1.0/(OMEGA^2)*(tabV[i] - tabV[j])*inv_norm
            Vdot_ij = - (tabX[i] - tabX[j])*inv_norm
            tabXdot[i] += Xdot_ij*tabm[j]
            tabVdot[i] += Vdot_ij*tabm[j]
            tabXdot[j] -= Xdot_ij*tabm[i]
            tabVdot[j] -= Vdot_ij*tabm[i]
        end

    end

end

"""
Function that compute the two body interraction for a given pair
# Arguments
+ index of the first particle
+ index of the second particle
+ dt time step making the rotation
"""
function pair_two_body!(index_first_particle::Int64, index_second_particle::Int64, dt::Float64)

    #index particules info at time 0
    X1 = tabX[index_first_particle]
    Y1 = tabV[index_first_particle]/OMEGA
    m1 = tabm[index_first_particle]
    X2 = tabX[index_second_particle]
    Y2 = tabV[index_second_particle]/OMEGA
    m2 = tabm[index_second_particle]

    #compute the reduce mass & the angular frequency of the pair interraction : O
    mu = (m1*m2)/(m1 + m2)
    O = (2*G*m1*m2)/(pi*mu*OMEGA)*(1.0/sqrt((X1 - X2)^2 + (Y1 - Y2)^2))

    #barycenter
    XG = (m1*X1 + m2*X2)/(m1 + m2)
    YG = (m1*Y1 + m2*Y2)/(m1 + m2)
    #complex barycenter
    ZG = XG + im*YG
    #initial complex position
    Z1 = X1 + im*Y1
    Z2 = X2 + im*Y2
    #initial deltaZ
    deltaZ = Z1 - Z2

    rot = exp(-im*O*dt)
    #turn of dt
    Z1 = ZG + (mu/m1)*deltaZ*rot
    Z2 = ZG - (mu/m2)*deltaZ*rot
    #update position
    tabX[index_first_particle] = real(Z1)
    tabV[index_first_particle] = imag(Z1)*OMEGA

    tabX[index_second_particle] = real(Z2)
    tabV[index_second_particle] = imag(Z2)*OMEGA
 
end

"""
Function that compute the two body interraction for a given pair
# Arguments
+ index of the first particle
+ index of the second particle
+ dt time step making the rotation
+ sign of the movment
"""
function pair_two_body_real!(index_first_particle::Int64, index_second_particle::Int64, dt::Float64)

    #index particles info at time 0
    X1 = tabX[index_first_particle]
    Y1 = tabV[index_first_particle]/OMEGA
    m1 = tabm[index_first_particle]
    X2 = tabX[index_second_particle]
    Y2 = tabV[index_second_particle]/OMEGA
    m2 = tabm[index_second_particle]

    #compute the reduce mass & the angular frequency of the pair interraction : O
    mu = (m1*m2)/(m1 + m2)
    O = ((2*G*m1*m2)/(pi*mu*OMEGA)*(1.0/sqrt((X1 - X2)^2 + (Y1 - Y2)^2)))

    #barycenter
    XG = (m1*X1 + m2*X2)/(m1 + m2)
    YG = (m1*Y1 + m2*Y2)/(m1 + m2)
   
    #initial deltaZ
    deltaX = X1 - X2
    deltaY = Y1 - Y2

    s, c = sincos(O*dt)

    X1 = XG + (mu/m1)*(deltaX*c + deltaY*s)
    Y1 = YG - (mu/m1)*(deltaX*s - deltaY*c)

    X2 = XG - (mu/m2)*(deltaX*c + deltaY*s)
    Y2 = YG + (mu/m2)*(deltaX*s - deltaY*c)

    #update position
    tabX[index_first_particle] = X1
    tabV[index_first_particle] = Y1*OMEGA

    tabX[index_second_particle] = X2
    tabV[index_second_particle] = Y2*OMEGA

end

"""
apply the backward interraction for dt
# Argument
+ dt, step time
"""
function backward!(dt::Float64)

    for n = 1:NPART

        X = tabX[n]
        Y = tabV[n]/OMEGA

        Z = X + im*Y

        Z *= exp(-im*(OMEGA/2)*dt)

        tabX[n] = real(Z)
        tabV[n] = imag(Z)*OMEGA

    end

end

"""
apply the backward interraction for dt
# Argument
+ dt, step time
"""
function backward_real!(dt::Float64)

    for n = 1:NPART

        X = tabX[n]
        Y = tabV[n]/OMEGA

        s, c = sincos((OMEGA/2)*dt)

        tabX[n] = X*c - Y*s
        tabV[n] = (X*s + Y*c)*OMEGA
        

    end

end

"""
Compute time derivative 
"""
function tab_dot!()

    if (MEAN)

        tab_dot_mean!()

    end

end