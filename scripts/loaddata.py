#!/usr/bin/env python
import numpy             as np
import matplotlib.pyplot as plt

import illstack as istk

import mpi4py.rc

istk.init.initialize('istk-params.txt')

field_list = ['Coordinates','Masses']
gas_particles = istk.io.getparticles(135,'gas',field_list)

field_list = ['GroupMass','GroupPos','Group_R_Mean200']
halos = istk.io.gethalos(135,field_list)

xp = gas_particles['Coordinates'][:,0]
yp = gas_particles['Coordinates'][:,1]
zp = gas_particles['Coordinates'][:,2]

xh = halos['GroupPos'][:,0]
yh = halos['GroupPos'][:,1]
zh = halos['GroupPos'][:,2]

rh = halos['Group_R_Mean200']
mh = halos['GroupMass']

if istk.params.rank==0:
    print 'min, max of halo     x-coordinates: ',xh.min(),xh.max()
    print 'min, max of particle x-coordinates: ',xp.min(),xp.max()

istk.cystack.cullonhalos(xp,yp,zp,xh,yh,zh,rh,mh)

