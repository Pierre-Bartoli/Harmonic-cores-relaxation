##################INSTANTANEOUS#################
"""Mass of particles
"""
const tabm      = zeros(Float64,NPART)

"""position of particles
"""
const tabx      = zeros(Float64,NPART)

"""velocities of particles
"""
const tabv      = zeros(Float64,NPART)

"""position in X
"""
const tabX      = zeros(Float64, NPART)

"""velocity in V
"""
const tabV      = zeros(Float64, NPART)

"""acceleration of particles
"""
const tabvdot   = zeros(Float64,NPART)

################CUMULATIVE########################

##################################################
# Table of the respective order of the particles
# The convention is that tabperm[i] corresponds
# to the index of the particle which is at the i-th
# position starting from the left
"""Table of the particles' order from left to right
"""
const tabperm = zeros(Int64, NPART)

#####
# We initialise the permutation array
# so that is already a permutation of 1:N
# from the start
for i = 1:NPART
    tabperm[i] = i
end

##################################################
# Table of the cumulative masses from left to right
# The mass of the particle at i-th postition
# is included in tabcummm[i]
# As a consequence, tabcumm[NPART] = MTOT

"""Table of the cumulative mass from left to right
"""
const tabcumm = zeros(Float64, NPART)

##################################################
# Table of the cumulative masses from right to left
# The mass of the particle at i-th postition
# is included in tabcummm[i]
# As a consequence, tabcumm[1] = MTOT

"""Table of the cumulative mass from right to left
"""
const Q = zeros(Float64, NPART)
"""Table of the cumulative mass time position from right to left
"""
const P = zeros(Float64, NPART)

"""
This function fill the array tabX and tabV 
# Argument 
+ time when tabX and tabV is needed
"""
function fill_tabXV!(time::Float64)

    #fill tabX and tabV
    for n=1:NPART

        X, V = xv_to_XV(tabx[n], tabv[n], time)

        tabX[n] = X
        tabV[n] = V

    end

end
"""
Computes the total momentum of a system
Ptot = SUM[m_i*v_i,{i,1,NPART}]
"""
function get_P_tot()
    Ptot = 0.0 # Initialising the counter
    #####
    for i=1:NPART # Loop over the particles
        Ptot += tabm[i]*tabv[i] # Adding the momentum of the particle i
    end
    #####
    return Ptot # Output
end

"""
Computes the total energy of a system \n
Etot = E_grav + E_kin \n
E_kin = SUM[(1/2)m_i*v_i^2,{i,1,Npart}] \n
E_grav = SUM[m_i*m_j*G*|x_i-x_j|,{i<j}]
"""
function get_E_tot()

    #sort the particles
    tabperm!()

    q = 0.0
    p = 0.0
    #compute sum in reverse (indices are N - n)
    for n=1:NPART # it is sum_j=i mj
        ind = tabperm[NPART - n + 1] # get the indices at the good emplacement
        q += tabm[ind] # add cumulative mass
        p += tabm[ind]*tabx[ind] # cumulative tabm*tabx
        Q[NPART - n + 1] = q # put in array in reverse side
        P[NPART - n + 1] = p # put in array in reverse side
    end

    E_grav = 0.0
    for n=1:NPART #then compute the energy
        ind = tabperm[n] # get the right index
        S = tabm[ind]*(tabx[ind]*Q[n] - P[n])# as in notes
        E_grav -= G*S #we are going in reverse wise
    end

    E_kin = 0.0 # energy kin
    for i=1:NBATHPART
        E_kin += 0.5*tabv[i]^2*tabm[i] # E_kin = 1/2mv^2
    end
    
    return E_grav + E_kin
end