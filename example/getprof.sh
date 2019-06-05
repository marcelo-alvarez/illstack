#!/bin/bash
#Columns: 1 init, 2 profile type, 3 lower mass limit, 4 upper mass limit, 5 snapshot (redshift) 
#python ../scripts/profiles.py istk-params.txt dmdens 8e12 1.2e13 084 # dm
#python ../scripts/profiles.py istk-params.txt gasdens 8e11 1.2e12 120 # gas
#python ../scripts/profiles.py istk-params.txt gaspth 8e11 1.2e12 120   # gas
python ../scripts/showstacks_profiles.py stack_gasdens_tng_12.npz stack_gaspth_tng_12.npz 8e11 1.2e12 dens_av
