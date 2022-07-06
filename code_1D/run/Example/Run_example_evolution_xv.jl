##################################################
include("../../src/Main.jl") # All the include
##################################################

################CONST#############################

"""x_t in time of the simulation
"""
const x_t = zeros(Float64, NPART, NDUMP + 1)

"""v_t in time of the simulation
"""
const v_t = zeros(Float64, NPART, NDUMP + 1)

"""Time @ istep
"""
const t = zeros(Float64, NDUMP + 1)

##################################################

"""
Dump the needed information in the arrays
"""
function dump!(istep)

    #fill tabx and tabv
    for n=1:NPART
        x_t[n, istep] = tabx[n]
        v_t[n, istep] = tabv[n]
    end

    #fill the time
    t[istep] = get_time()
end

"""
Launch the simulation, integrate x and v during the motion
"""
function run!()

    dump!(1)#dump a time 0

    for idump=1:NDUMP

        integrate_NSTEPS!() # Integrating for a few timesteps
        dump!(idump + 1)# dump at every dump

    end
end

"""
Function to perform a dump, it write all the data of the simulation in arrays of size NDUMP
# Argument
nameFile the name of the file that dump the simulation
"""
function write!(namefile::String)

    # suppress the previous file
    if(isfile(namefile))# but to supress it be sure that it exist
        rm(namefile) #remove it
    end

    file = h5open(namefile,"w") # Opening the file

    #########################VARIABLES####################

    write(file, "x_t", x_t)
    write(file, "v_t", v_t)
    write(file, "t", t)

    #########################ARGS##########################

    write(file, "NTESTPART", NTESTPART)
    write(file, "NBATHPART", NBATHPART)
    write(file, "DT", DT)
    write(file, "NDUMP", NDUMP)
    write(file, "NSTEPS", NSTEPS)
    write(file, "A", A)
    write(file, "SELF", SELF)
    write(file, "HARM", HARM)
    write(file, "SEED", SEED)

    ########################CONST###########################

    write(file, "NPART", NPART)
    write(file, "MTOT", MTOT)
    write(file, "OMEGA", OMEGA)
    write(file, "G", G)

    #########################################################

    close(file) # Closing the file
    
end

"""
Give an output of the simulation (eg. video, plot...)
"""
function output!(path)

    #output a mp4 
    anim = @animate for i=1:NDUMP
        #scatter particles for the animation
        scatter(x_t[:, i], 
            v_t[:, i],
            xlabel = "x", 
            ylabel = "v", 
            label = "", 
            size = (500, 500), 
            title = "Evolution of x, v for \n t = "*string(round(t[i]*(OMEGA/2pi), digits=2))*" Td",
            xlim = (-1.5, 1.5),
            ylim = (-1.5, 1.5),
            aspect_ratio = :equal)
    end

    gif(anim, path, fps = 30)

end
##################################################
init_sampling!() # Sampling the initial location of the particles
##################################################
println("Npart | ", NPART, " | DT | ", DT, " | Nsteps | ", NSTEPS, "| NDUMP| ", NDUMP, " | seed | ", SEED)
##################################################
print("running : ")
@time run!() #run the simulation
println("Numerical simulation finished")
##################################################
if (WRITE)
    print("write : ")
    @time write!(PATH*".hf5") # Dumping
    println("Writing finished")
end
##################################################
print("output : ")
@time output!(PATH*".mp4")
println("Succesfully finished")

