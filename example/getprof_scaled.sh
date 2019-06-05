#!/bin/bash
#Columns: 1 init, 2 profile type, 3 lower mass limit, 4 upper mass limit, 5 snapshot (redshift) 
python ../scripts/profiles_scaled.py istk-params_scaled.txt dmdens 8e11 1.2e12 084 # dm
python ../scripts/profiles_scaled.py istk-params_scaled.txt gasdens 8e11 1.2e12 084 # gas
python ../scripts/profiles_scaled.py istk-params_scaled.txt gaspth 8e11 1.2e12 084   # gas
python ../scripts/showstacks_scaled.py stack_gasdens_scaled_12.npz stack_gaspth_scaled_12.npz 8e11 1.2e12 dens
