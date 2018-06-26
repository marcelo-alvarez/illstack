#!/usr/bin/env python
import sys

import numpy             as np
import matplotlib.pyplot as plt
import illstack as istk

import mpi4py.rc

istk.init.initialize(sys.argv[1])

# here we want gas density so the value to bin is mass with a volume weight = true

field_list = ['Coordinates','Masses']
gas_particles = istk.io.getparticles(135,'gas',field_list)

field_list = ['GroupMass','GroupPos','Group_R_Mean200']
if istk.params.rank==0: print 'reading halos'
halos = istk.io.gethalos(135,field_list)

posp = gas_particles['Coordinates'] 
vals = gas_particles['Masses'] # here we want gas density so the value to bin is mass with a volume weight = true

posh = halos['GroupPos']
mh   = halos['GroupMass']
rh   = halos['Group_R_Mean200']

ntile = 10
volweight = True # here we want gas density so the value to bin is mass with a volume weight = true
r, val, n, mh = istk.cystack.stackonhalos(posp,vals,posh,mh,rh,
                                         ntile,volweight)
print 'shapes: ',np.shape(r),np.shape(val),np.shape(n),np.shape(mh)

np.savez('stack.npz',r=r,val=val,n=n,mh=mh)
