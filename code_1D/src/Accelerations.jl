"""
Sorts the particles according to their position \n
The convention is that tabperm[i]
returns the index of the particle that is at the i-th
position starting from the left
"""
function tabperm!()
    # We sort the particles according to their position
    # The syntax used allows for a minimal number of allocations
    sort!(tabperm,QuickSort,Perm(ForwardOrdering(),tabx))
end

"""
Computing the cumulative mass \n
The convention is that tabcumm[i]
is the mass that is on the left of the i-th particle
included its own individual mass \n
As a consequence, tabcumm[NPART] = Mtot
"""
function tabcumm!()
    # First, we deal with the first particle on the left
    ind_loc = tabperm[1] # Physical index of the particle that is the most on the left
    cumm_loc = tabm[ind_loc] # Initialising the cumulative mass with the individual mass of the left-most particle
    #####
    tabcumm[1] = cumm_loc # Filling in the array of the cumulative mass
    #####
    for i=2:NPART # Loop over the rest of the particles. ATTENTION, the loop starts at 2
        ind_loc = tabperm[i] # Physical index of the particle that is at the i-th position from the left
        cumm_loc += tabm[ind_loc] # Adding the mass of the i-th particle from the left
        #####
        tabcumm[i] = cumm_loc # Filling in the array of the cumulative mass
        #####
    end
end

##################################################
# Computes the accelerations of the particles
# The acceleration, vdot, on a particle is
# vdot = G*(M_right - M_left)
#      = G*(Mtot - 2*M_left - m)  [using M_right + M_left = MTOT - m]
#      = G*(Mtot - 2*Mcum + m)    [using M_left = M_cum - m]
##################################################

"""
Computes the accelerations of the particles \n
The acceleration, vdot, on a particle is \n
vdot = G*(M_right - M_left) \n
= G*(Mtot - 2*M_left - m)  [using M_right + M_left = MTOT - m] \n
= G*(Mtot - 2*Mcum + m)    [using M_left = M_cum - m] \n
"""
function tabvdot_self!()
    #####
    tabperm!() # We sort the particles
    #####
    tabcumm!() # We compute the cumulative mass
    #####
    for i=1:NPART # Loop over the particles from left to right
        #####
        ind_loc = tabperm[i] # Physical index of the particle that is at the i-th position from the left
        #####
        cumm_loc = tabcumm[i] # Cumulative mass on the left of the i-th particle
        m_loc = tabm[ind_loc] # Individual mass of the particle
        #####
        vdot = G*(MTOT - 2.0*cumm_loc + m_loc)# Acceleration of the i-th particle from the left
        #####
        tabvdot[ind_loc] += vdot # Updating the velocity of the particle
        #####
    end
end

"""
Add a harmonic force to the particle
"""
function tabvdot_harm!()
    for i=1:NPART
        tabvdot[i] -= OMEGA^2*tabx[i]
    end
end

"""
Compute the vdot of particles at given time
"""
function tabvdot!()
    #####
    fill!(tabvdot, 0.0) # Initialising the velocities 
    #####
    if(SELF)
        tabvdot_self!() # Compute the self gravity on the particles
    end
    #####
    if(HARM)
        tabvdot_harm!()
    end
end
