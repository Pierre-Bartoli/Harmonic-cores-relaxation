"""
Recenter the sampling, du to the random distribution all x and v are not centered anymore
"""
function center!()

    Xg_loc = 0.0
    Vg_loc = 0.0
    
    for n=1:NPART
        Xg_loc += tabm[n]*tabX[n]
        Vg_loc += tabm[n]*tabV[n]
    end
    Xg_loc /= MTOT
    Vg_loc /= MTOT

    for n=1:NBATHPART
        tabX[n] -= Xg_loc
        tabV[n] -= Vg_loc
    end
end

"""
Get the sample of bath paticles according to a constant density
# return
couple (X, V) sample to be in quasi stationnary state
"""
function getSample!()

    X = 2.0*A*rand() - A # Drawing x uniformly in [-1,1]
    V = -sqrt(G*MTOT*A - (OMEGA^2)*(X^2))*cos(rand()*pi) # Drawing v
    return [X, V]

end

"""
fill the particle of the bath
# Argument
+ n the n-th particle \n
+ m the mass of the particle \n
+ [X, V] the couple of the particle position
"""
function fillBathParticles!(n::Int64, m::Float64, coord::Array{Float64, 1})
    tabm[n] = m # Filling in the mass
    tabX[n] = coord[1] # Filling in the position
    tabV[n] = coord[2] # Filling in the momentum = m*v
end

"""
fill with test particles in a small patch
# Argument
+ n index of the tab
+ r radius of the patch
+ xC, yC center of the patch
"""
function fillTestParticles!(n, r, xC, yC)
    tabm[n] = 0.0 # filling with null mass
    X = 2.0*r*rand() - r + xC
    tabX[n] = X # filling uniformly at given point
    Y = 2.0*r*rand() - r + yC
    tabV[n] = Y # filling uniformly at given point
end

"""
Draws the initial statistics of the particles, to sample the harmonic distribution
"""
function init_sampling!()

    for n=1:(NBATHPART) # Loop over the particles
        m = MTOT/NBATHPART # Individual mass of the particles
        fillBathParticles!(n, m, getSample!()) # fill the bath particles
    end

    for n=NBATHPART+1:NPART # test particles
        fillTestParticles!(n, R, XC, YC) # fill within a circle of radius 0.01 and of center (0.1, 0.5)
    end

    center!() # Recenter the particles

end
