##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin

    "--path"
    help = "path of the saved file"
    arg_type = String
    default = "data/default/default"
    "--write"
    help = "if the program write a dump"
    arg_type = Bool
    default = false


    "--Nbathpart"
    help = "Number of particles"
    arg_type = Int64
    default = 20000
    "--xC"
    help = "center of test particles in coord X"
    arg_type = Float64
    default = 0.0
    "--yC"
    help = "center of test particles in coord Y"
    arg_type = Float64
    default = 0.0
    "--r"
    help = "radius of test particles"
    arg_type = Float64
    default = 0.0
    "--Ntestpart"
    help = "Number of particles"
    arg_type = Int64
    default = 1

    "--dt"
    help = "Integration timestep"
    arg_type = Float64
    default = 0.001
    "--Nsteps"
    help = "Number of steps dt used in integrate_Nsteps!()"
    arg_type = Int64
    default = 10
    "--Ndump"
    help = "Number of steps between dumps"
    arg_type = Int64
    default = 100
    "--integration"
    help = "integration method"
    arg_type = String
    default = "RK2"

    "--Nj"
    help = "Number of bin in action direction"
    arg_type = Int64
    default = 101
    "--Ntheta"
    help = "Number of bin in angle direction"
    arg_type = Int64
    default = 101

    "--self"
    help = "if the self gravity must be turned on"
    arg_type = Bool
    default = true
    "--harm"
    help = "if the harmonic force must be turned on"
    arg_type = Bool
    default = false
    "--mean"
    help = "if the Gamma is mean"
    arg_type = Bool
    default = true

    "--a"
    help = "limit of core (density function)"
    arg_type = Float64
    default = 1.0
    
    "--Nseed"
    help = "number of different seed for the random generator"
    arg_type = Int
    default = 1000
    "--seed"
    help = "Seed for the random generator"
    arg_type = Int
    default = 1
end
##################################################
parsed_args = parse_args(tabargs)
##################################################
# General parameters
##################################################

#######################PATH#######################
"""Path of the data
"""
const PATH = parsed_args["path"]

"""if the program write the data out
"""
const WRITE = parsed_args["write"]

#######################PARTICLES##################

"""Number of test particles
"""
const NTESTPART = parsed_args["Ntestpart"]

"""center of test particles in coord X
"""
const XC = parsed_args["xC"]

"""center of test particles in coord Y
"""
const YC = parsed_args["yC"]

"""radius of test particles in coord X
"""
const R = parsed_args["r"]

"""Number of bath particles
"""
const NBATHPART = parsed_args["Nbathpart"]

"""Number of particles
"""
const NPART = NTESTPART + NBATHPART

#####################TIME#########################

"""Integration timestep
"""
const DT = parsed_args["dt"]

"""Number of step between dumps
"""
const NDUMP = parsed_args["Ndump"]

"""Number of integration time steps
"""
const NSTEPS = parsed_args["Nsteps"]

"""How the movment is integrate
"""
const INTEGRATION = parsed_args["integration"]

#######################SIZE########################

"""Size of the core
"""
const A = parsed_args["a"]

######################BIN###########################

"""Number of bin in action
"""
const NJ = parsed_args["Nj"]

"""Number of bin in angle
"""
const NTHETA = parsed_args["Ntheta"]

"""Sample the bin in action
"""
const HJ = (A^2*0.5)/(NJ - 1)

"""Sample the bin in angle
"""
const HTHETA = 2pi/(NTHETA - 1)

#######################FORCES#######################

"""If the self gravity is turned on
"""
const SELF = parsed_args["self"]

"""If the harmonic force is turned on
"""
const HARM = parsed_args["harm"]

"""If use Gamma mean
"""
const MEAN = parsed_args["mean"]

#####################OTHERS#########################

"""Seed of the random
"""
const SEED = parsed_args["seed"]

"""Number of seed use
"""
const NSEED = parsed_args["Nseed"]

##################################################
Random.seed!(SEED) # Setting the random seed
##################################################
