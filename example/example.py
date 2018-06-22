#!/usr/bin/env python
import numpy             as np
import matplotlib.pyplot as plt

import illstack as istk
import illstack.snapshotinfo as ilsp

import mpi4py.rc

istk.init.initialize('istk-params_template.txt')

mean_gas_mass = ilsp.meanparticleproperty(135,'gas','Masses')
print('mean gas particle mass = ',mean_gas_mass)

