#!/bin/bash
PATHFILE=data/example/example_evolution_xv
WRITE=false

NBATHPART=2000 #Number of bath particles
NTESTPART=0 #Number of test particles

DT=0.001 #10^-3 TD
NSTEPS=50
NDUMP=100

SELF=true #only self gravitation
HARM=false #no harmonic force

A=1 #size of the core
SEED=1 #seed of the realisation

JULIA=julia
RUN=run/Example/Run_example_evolution_xv.jl
PREFIX= #directory of the code
LOG=log/log_example_evolution_xv

cd ${PREFIX}
${JULIA} ${RUN} --path ${PATHFILE} --write ${WRITE} --Nbathpart ${NBATHPART} \
--Ntestpart ${NTESTPART} --dt ${DT} --Nsteps ${NSTEPS} \
--Ndump ${NDUMP} --self ${SELF} --harm ${HARM} \
--a ${A} --seed ${SEED} | tee ${LOG}