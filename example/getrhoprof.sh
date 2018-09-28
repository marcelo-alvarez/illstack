#!/bin/bash

python ../scripts/masscutrho.py istk-params.txt 0 8e12     # dm
python ../scripts/masscutrho.py istk-params.txt 1 8e12     # gas
python ../scripts/showstacks.py stack_dm.npz stack_gas.npz # show the stacks
