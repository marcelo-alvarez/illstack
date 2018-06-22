#!/usr/bin/env python
import numpy             as np
import matplotlib.pyplot as plt

import illstack as istk
import illstack.snapshotinfo as ilsp

import mpi4py.rc

istk.init.initialize('istk-params.txt')

field_list = ['Coordinates','Masses']
gas_particles = istk.io.getparticles(135,'gas',field_list)

field_list = ['GroupMass','GroupPos','Group_R_Mean200']
halos = istk.io.gethalos(135,field_list)

xp = gas_particles['Coordinates'][:,0]
yp = gas_particles['Coordinates'][:,1]
zp = gas_particles['Coordinates'][:,2]

print np.shape(xp)


