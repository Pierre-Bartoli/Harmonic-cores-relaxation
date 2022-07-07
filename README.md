# Harmonic-cores-relaxation
This code computes the motion of N self-gravitating particles in a 1D harmonic core (i.e. a constant density profile). For constant density profile in 1D, all the particles in the **physcial space** have the same frequency for their orbit (time to make an orbit is the same for every particle). In **phase space**, this motion corresponds to an harmonic oscillator (particles are turning around the center of phase). On longer times, these circles are slowly deformed. To better track these deformations, we remove the *average harmonic motion* by performing a change of variables, **(x, v)** -> **(X, V)**, where x stands for position and v for the velocity of the particle. 

There are two codes for this this N-body problem. First the **code_1D** computes the dynamics directly in **(x, v)**. Second, **code_1D_XV** computes the dynamics within the **(X, V)** variables, with a time-averaged interaction potential. We note that **code_1D** propagates ONLY **(x, v)** but can render the **(X, V)** motion by doing a change of variables afterwards. HOWEVER this dynamics is NOT averaged on time (that's why particles "shaking" in **code_1D**). Inversly, **code_1D_XV** computes ONLY the **(X, V)** motion AND is averaged on time.

What we can get:

![evolution_xv](https://user-images.githubusercontent.com/108795620/177579142-d092ceec-4d4e-4111-b980-bbd4a26481ca.png)


## Code_1D
**Code_1D** integrates the dynamics using a classical leapfrog integrator, which is symplectic. The code operates on the dynamical time TD (orbital time) and typically, one picks 10<sup>-3</sup> TD as the time step DT. The complexity of the algorithm is O(NlnN).

## Code_1D_XV
**Code_1D_XV** integrates the dynamics with different integrators, such as Runge-Kutta or Gauss-Legendre. The code operates on the balistic time TB (time of deformation, TB = sqrt(N)TD). Hence timesteps are larger than TD, BUT the integrator complexity is O(N<sup>2</sup>).

# How to use the code
The code is in julia. To run it, you must install first some packages.

## Packages
To install packages, launch julia and type:

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

![example_evolution_xv](https://user-images.githubusercontent.com/108795620/177580140-7b0c4991-2be8-42e1-be9d-0b659f0b5e74.gif)


- To run the provided examples, download the code, go to `job/Example` folder.

- Then give access to execute: `chmod +x job_example_evolution_xv.sh`

- Finally launch the bash file: `.\job_example_evolution_xv.sh`, the output will be in *data/example* and a log file will be located in *log/example*.

## Output others parameters

The preferable way to do it is: 
- Write your *Run.jl* in Run folder
- Write *job.sh* in job folder to output the wanted parameters.
- Do not forget to change **RUN**, **LOG** and **PATHFILE** argument in *job.sh* (where path file is the output file).
