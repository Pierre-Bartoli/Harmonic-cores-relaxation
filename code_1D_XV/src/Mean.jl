"""Acceleration of the X"""
const tabXdot = zeros(Float64, NPART)
"""Acceleration of the V"""
const tabVdot = zeros(Float64, NPART)

"""Inverse of the norm to compute this amount once"""
const inv_norm = zeros(Float64, NPART)

"""V particles positions"""
const tabV = zeros(Float64, NPART)
"""X particles positions"""
const tabX = zeros(Float64, NPART)
"""mass of particles"""
const tabm = zeros(Float64, NPART)

"""
This function compute the total energy of the bath
# return
+ Energy of the bath (needs to be invariant)
"""
function E_tot()

    E_kin = 0

    #compute the backward
    for i = 1:NPART
        E_kin -= 0.25*(tabm[i]*tabV[i]^2 + tabm[i]*tabX[i]^2*OMEGA^2)
    end

    E_grav = 0

    #compute the coupled
    for i = 1:NPART
        for j = (i+1):NPART

            E_grav += tabm[i]*tabm[j]*(2.0*G/(pi))*sqrt((tabX[i] - tabX[j])^2 + (1/OMEGA)^2*(tabV[i] - tabV[j])^2)

        end
    end

    return E_kin + E_grav

end

function B_tot(t::Float64)

    A = 0

    for n=1:NPART

        A += tabm[n]*tabX[n] + (im*tabm[n]*tabV[n])/OMEGA

    end

    B_tot = A*exp(im*0.5*OMEGA*t)

    return abs(B_tot)

end

function I_tot()

    I = 0

    for n=1:NPART
        I += tabm[n]^2*OMEGA^2*tabX[n]^2 + tabm[n]^2*tabV[n]^2
    end

    return I

end
