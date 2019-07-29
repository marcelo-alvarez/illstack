#!/bin/bash
#Columns: 1 init, 2 profile type, 3 lower mass limit, 4 upper mass limit, 5 snapshot (redshift) 
#python ../scripts/profiles_colors.py istk-params.txt dmdens 8e12 1.2e13 084 # dm
#python ../scripts/profiles_colors.py istk-params.txt gasdens 8e12 1.2e13 084 # gas
#python ../scripts/profiles_colors.py istk-params.txt gaspth 8e12 1.2e13 084   # gas
python ../scripts/showstacks_colors.py stack_gasdens_red_12.npz stack_gaspth_red_12.npz stack_gasdens_blue_12.npz stack_gaspth_blue_12.npz 8e11 1.2e12 pres
