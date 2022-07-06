#!/bin/bash
PATHFILE=data/example/example_evolution_XV
WRITE=false

INTEGRATION=GL2 #integration scheme, here runge-kutta order 4

NBATHPART=2000 #Number of bath particles
NTESTPART=0 #Number of test particles

DT=1 #10^1 TD
NSTEPS=2
NDUMP=100

SELF=true #only self gravitation
HARM=false #no harmonic force

A=1 #size of the core
SEED=1 #seed of the realisation

JULIA=julia
RUN=run/Example/Run_example_evolution_XV.jl
PREFIX= #directory of the code
LOG=log/log_example_evolution_XV

cd ${PREFIX}
${JULIA} ${RUN} --path ${PATHFILE} --write ${WRITE} --Nbathpart ${NBATHPART} \
--Ntestpart ${NTESTPART} --dt ${DT} --Nsteps ${NSTEPS} \
--Ndump ${NDUMP} --self ${SELF} --harm ${HARM} \
--a ${A} --seed ${SEED} --integration ${INTEGRATION} | tee ${LOG}
