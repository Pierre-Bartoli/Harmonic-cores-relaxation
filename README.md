# Harmonic-cores-relaxation
This code compute the motion of N particles in an harmonic cores (i.e. constant density profile) in 1D. All particles are orbiting and the interaction is the gravity of every other particles. For constant density profile in 1D, we can show that every particles in the **physcial space** have the same frequency for their orbit (time to make an orbit is the same for every particles). In the **phase space** this motion correspond to an oscillator (particles are turning around the center) but with lots of time those circles in the phase space are deformed. To have a better view of the deformation, one can remove the *average harmonic motion* by performing a variable change, **(x, v)** -> **(X, V)**, where x stand for position and v for the velocity of the particle. 

There is two code to compute this N-body problem. First the **code_1D** compute the real problem directly in **(x, v)** the true variable that appear in "real life". And then the **code_1D_XV** compute the **(X, V)** motion, without the *average harmonic motion* and with a time average. Notice that the **code_1D** propagate ONLY the **(x, v)** but can render the **(X, V)** motion by doing a variable change after having computed the **(x, v)** motion HOWEVER the motion is NOT average on time (that why *example* in phase space are "shaking" in **code_1D**). Inversly **code_1D_XV** compute ONLY the **(X, V)** motion AND is  average on time, so it cannot go back the true **(x, v)**.

## Code_1D
**Code_1D** propagate the motion with a classical leapfrog integrator, which is sympletic. The code must operate on dynamical time TD (time to do an orbit) and typically to avoid error, one must put 10<sup>-3</sup> TD as time step DT. The complexity of the algorithm is O(NlnN).

## Code_1D_XV
**Code_1D_XV** propagate the motion with different integrator as Runge-Kutta or Gauss-Legendre. The code operates on balistic time TB (time of deformation, TB = sqrt(N)TD) which is larger than TD, so integration can be faster BUT the integrator complexity is O(N<sup>2</sup>).

# How to use the code
The code is in Julia, to make it run, you must install first some packages.

## Packages
To install package you have to launch Julia, and type in :

```
using Pkg
Pkg.add("Package Name")
```

You must install :
- ArgParse, to parse command-line arguments
- BenchmarkTools, to have good timing measurements
- Random, to be able to draw random numbers
- HDF5, to use HDF5 files
- Plots, to use plots
by doing,
```
using Pkg
Pkg.add("ArgParse")
Pkg.add("BenchmarkTools")
Pkg.add("Random")
Pkg.add("HDF5")
Pkg.add("Plots")
```

## Run examples

- Create a directory name log in **code_1D** and in **code_1D_XV** folder.

- To run examples you must download the code, go to job folder and open one of *job_example_evolution*. Then you must remplace the comment at **PREFIX** line to put your code_1D (or code_1D_XV) path. 

For example, if you want to run the *job_example_evolution_xv.sh*, you must find the line **20** and put 
`PREFIX=home/user/Download/code_1D/` if the code_1D is located in Download.

- Then launch the bash file : `.\job_example_evolution_xv.sh`, the output will be in *data/example* and a log file will be located in *log/example*.

## Output others parameters

The preferable way to do it is : 
- you should write your *Run.jl* in Run folder
- and write *job.sh* in job folder to output the wanted parameters.
- don't forget to change **RUN**, **LOG** and **PATHFILE** argument in *job.sh* (where path file is the output file).

## Use written file (if WRITE value is true)

You may want write information on hf5 file, to do so it is easy, with **WRITE**=true in *job.sh*. You must write information needed in the Run.jl, in write! function:

```
function write!(namefile::String)

    # suppress the previous file
    if(isfile(namefile))# but to supress it be sure that it exist
        rm(namefile) #remove it
    end

    file = h5open(namefile,"w") # Opening the file

    write(file, "x_t", x_t) #put x_t array in "x_t"
    write(file, "v_t", v_t)
    write(file, "t", t)
    
    write(file, "value", value) #general value

    close(file) # Closing the file
    
end
```
And then you can read it with an other code. For example one can type in Julia this code :
```
using HDF5

fid = h5open("path_to_hf5_file.hf5","r")

#get x, v & time
x_t = read(fid["x_t"])
v_t = read(fid["v_t"])
t = read(fid["t"])

#general value (must be written in Run.jl)
value = read(fid["value"])
```
