#!/bin/bash

python ../scripts/masscutrho.py istk-params.txt 0 1e13     # dm
python ../scripts/masscutrho.py istk-params.txt 1 1e13     # gas
python ../scripts/showstacks.py stack_dm.npz stack_gas.npz # show the stacks
