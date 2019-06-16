#!/bin/bash
#Columns: 1 init, 2 profile type, 3 lower mass limit, 4 upper mass limit, 5 snapshot (redshift) 
#python ../scripts/profiles_scaled.py istk-params_scaled.txt dmdens 1e13 120 # dm
#python ../scripts/profiles_scaled.py istk-params_scaled.txt gasdens 1e13 120 # gas
#python ../scripts/profiles_scaled.py istk-params_scaled.txt gaspth 1e13 120   # gas
python ../scripts/showstacks_scaled.py stack_gasdens_scaled_ill_13.npz stack_gaspth_scaled_ill_13.npz 1e13 pres
