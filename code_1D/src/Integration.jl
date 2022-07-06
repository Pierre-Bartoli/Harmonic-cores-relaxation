"""Global time of the simulation
"""
const CURRENT_TIME = zeros(Int64, 1)

"""
Performs the integration of the particle's
using the standard leap-frog
"""
function integrate_DT!()
    #####
    # Step 1
    #####
    for n=1:NPART
        tabx[n] += 0.5*DT*tabv[n] # Drift for half a timestep
    end
    #####
    tabvdot!() # Computing the velocities
    #####
    for n=1:NPART
        tabv[n] += DT*tabvdot[n] # Kick for a full timestep + harmonic force
    end
    #####
    # Step 2
    #####
    for n=1:NPART
        tabx[n] += 0.5*DT*tabv[n] # Drift for half a timestep
    end
    #####

    CURRENT_TIME[1] += 1
end

"""
Wrapped function that performs the integration
for NSTEPS steps
"""
function integrate_NSTEPS!()

    for istep=1:(NSTEPS) # Loop over the steps to perform

        integrate_DT!() # Integrating for one timestep

    end

end

"""
Give the current time of the simulation
# Return 
+ current time
"""
function get_time()
    return CURRENT_TIME[1]*DT
end
