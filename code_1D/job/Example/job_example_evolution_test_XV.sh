#!/bin/bash
PATHFILE=data/example/example_evolution_test_XV
WRITE=false

NBATHPART=2000 #Number of bath particles
NTESTPART=100 #Number of test particles

XC=0.1 #position x of the blob
YC=0.1 #position y of the blob
R=0.05 #radius of the blob

DT=0.001 #10^-3 TD
NSTEPS=2000
NDUMP=100

SELF=true #only self gravitation
HARM=false #no harmonic force

A=1 #size of the core
SEED=1 #seed of the realisation

JULIA=julia
RUN=run/Example/Run_example_evolution_test_XV.jl
PREFIX= #directory of the code
LOG=log/log_example_evolution_test_XV

cd ${PREFIX}
${JULIA} ${RUN} --path ${PATHFILE} --write ${WRITE} --Nbathpart ${NBATHPART} \
--Ntestpart ${NTESTPART} --xC ${XC} --yC ${YC} --r ${R} --dt ${DT} \
--Nsteps ${NSTEPS} --Ndump ${NDUMP} --self ${SELF} --harm ${HARM} \
--a ${A} --seed ${SEED} | tee ${LOG}