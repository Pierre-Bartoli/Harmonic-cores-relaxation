# Physical constants

"""Gravitational constant. The pairwise interaction is U(x,x')=G|x-x'|
"""
const G = 1.0

"""Total mass of the system
"""
const MTOT = 1.0

"""Frequency of harmonic oscillator
"""
const OMEGA = sqrt((G*MTOT)/A)

"""Energy zero of the system
"""
const E_0 = G*MTOT*A