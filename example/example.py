#!/usr/bin/env python
import numpy             as np
import matplotlib.pyplot as plt

import illstack as istk
#import illstack.snapshotinfo as ilsp
from illstack.CompHaloProperties import CompHaloProp

import mpi4py.rc

istk.init.initialize('istk-params.txt')

#mean_gas_mass = ilsp.meanparticleproperty(135,'gas','Masses')
#print('mean gas particle mass = ',mean_gas_mass)


lims = [0.1,10.]
bins = 10
CHP = CompHaloProp(lims,bins)


pos = np.exp(np.random.rand(3,100)*2.0)

print pos

temp = (np.random.rand(100))
weight = np.array(1.+temp*0.0)


print temp.shape, np.mean(weight)

print CHP.radbins, CHP.BinCenter
print CHP.ComputeHaloProfile(pos,temp,weight)
