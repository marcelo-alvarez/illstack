#!/bin/bash

python ../scripts/profiles.py istk-params.txt dmdens 1e13     # dm
python ../scripts/profiles.py istk-params.txt gasdens 1e13     # gas
python ../scripts/profiles.py istk-params.txt gaspth 1e13     # gas
python ../scripts/showstacks.py stack_gasdens.npz stack_gaspth.npz # show the stacks
