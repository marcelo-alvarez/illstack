#!/bin/bash
#Columns: 1 init, 2 profile type, 3 lower mass limit, 4 upper mass limit, 5 snapshot (redshift) 
python ../scripts/profiles.py istk-params.txt dmdens 8e12 1.2e13 120 # dm
python ../scripts/profiles.py istk-params.txt gasdens 8e12 1.2e13 120 # gas
python ../scripts/profiles.py istk-params.txt gaspth 8e12 1.2e13 120   # gas
python ../scripts/showstacks_mean_profiles.py stack_gasdens.npz stack_gaspth.npz 8e12 1.2e13
