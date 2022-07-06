"""
Recenter the sampling, du to the random distribution all x and v are not centered anymore
"""
function center!()

    #init xg & vg
    xg_loc = 0.0
    vg_loc = 0.0
    
    #compute xg & vg
    for n=1:NPART
        xg_loc += tabm[n]*tabx[n]
        vg_loc += tabm[n]*tabv[n]
    end
    xg_loc /= MTOT
    vg_loc /= MTOT

    #apply it to the coords
    for n=1:NBATHPART
        tabx[n] -= xg_loc
        tabv[n] -= vg_loc
    end
end

"""
Get the sample of bath paticles according to a constant density
# return
couple (x, v) sample to be in quasi stationnary state
"""
function getSample()

    x = 2.0*A*rand() - A # Drawing x uniformly in [-1,1]
    v = -sqrt(G*MTOT*A - (OMEGA^2)*(x^2))*cos(rand()*pi) # Drawing v
    return [x, v]

end

"""
fill the particle of the bath
# Argument
+ n the n-th particle \n
+ m the mass of the particle \n
+ [x, v] the couple of the particle position
"""
function fillBathParticles!(n::Int64, m::Float64, coord::Array{Float64, 1})
    tabm[n] = m # Filling in the mass
    tabx[n] = coord[1] # Filling in the position
    tabv[n] = coord[2] # Filling in the velocity
end

"""
fill with test particles in a little patch
# Argument
+ n index of the tab
+ r radius of the patch
+ xC, yC center of the patch
"""
function fillTestParticles!(n, r, xC, yC)
    tabm[n] = 0.0 # filling with null mass
    x = 2.0*r*rand() - r + xC
    tabx[n] = x # filling uniformly at given point
    y = 2.0*r*rand() - r + yC
    tabv[n] = y # filling uniformly at given point
end

"""
Draws the initial statistics of the particles, to sample the harmonic distribution
"""
function init_sampling!()

    for n=1:(NBATHPART) # Loop over the particles
        m = MTOT/NBATHPART # Individual mass of the particles
        fillBathParticles!(n, m, getSample()) # fill the bath particles
    end

    for n=NBATHPART+1:NPART # test particles
        fillTestParticles!(n, R, XC, YC) # fill within a circle of radius 0.01 and of center (0.1, 0.5)
    end

    center!() # Recenter the particles

end
