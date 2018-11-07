#!/usr/bin/env python
import sys

import numpy             as np
import matplotlib.pyplot as plt
import illstack as istk

import mpi4py.rc

istk.init.initialize(sys.argv[1])
prof = str(sys.argv[2])  #This is getting specific profile (gas density, gas thermal pressure, dm)
mlow   = float(sys.argv[3])
mhigh = float(sys.argv[4])
ntile = 3 # controls tiling -- optimal value not yet clear
mhmin = mlow / 1e10 # minimum mass in 1e10 Msun/h
mhmax = mhigh / 1e10 # maximum mass in 1e10 Msun/h
volweight = True # here we want density so value to bin is mass with a volume weight = true
snap_num= float(sys.argv[5])  ##New, not hardcoded for z=0

omegadm = 0.3 - 0.042
#Xh=0.76
gamma = 5./3.

if prof=='gasdens':
    part_type='gas'
    field_list = ['Coordinates','Masses']
    gas_particles = istk.io.getparticles(snap_num,part_type,field_list) #Change redshift
    posp = gas_particles['Coordinates'] #position 
    vals = gas_particles['Masses'] 
elif prof=='dmdens':
    part_type='dm'
    # HARD CODED BOX SIZE 7.5e4 kpc/h
    part_massf=2.775e2*omegadm*(7.5e4)**3/1e10 # particle mass in 1e10 Msun/h
    field_list = ['Coordinates']
    posp = istk.io.getparticles(snap_num,part_type,field_list) #Change redshift
    vals = posp[:,0]*0 + part_massf / np.shape(posp)[0]
    print 'setting dm particle mass to = ',vals[0]*1e10,'Msun/h'
elif prof=='gaspth':
    part_type='gas'
    field_list = ['Coordinates','Masses','InternalEnergy']
    gas_particles = istk.io.getparticles(snap_num,part_type,field_list) #Change redshift
    posp = gas_particles['Coordinates']
    vals = gas_particles['Masses']*gas_particles['InternalEnergy']*(gamma-1.)
else:
    print "Please enter an appropriate option for the profile"
    print "gasdens,dmdens,gaspth"


field_list = ['GroupMass','GroupPos','Group_R_Mean200']
halos = istk.io.gethalos(snap_num,field_list) #Change redshift

posh = halos['GroupPos']
mh   = halos['GroupMass']
rh   = halos['Group_R_Mean200']

r, val, n, mh, rh, nprofs = istk.cyprof.stackonhalos(posp,vals,posh,mh,rh,
                                                  ntile,volweight,mhmin)
r  =np.reshape(r,  (nprofs,istk.params.bins))
val=np.reshape(val,(nprofs,istk.params.bins)) 
n  =np.reshape(n,  (nprofs,istk.params.bins))

print 'shapes: ','r', np.shape(r),'val', np.shape(val),'n', np.shape(n),'mh', np.shape(mh)

np.savez('stack_'+prof+'.npz',r=r[0],val=val,n=n,mh=mh,rh=rh,nprofs=nprofs,nbins=istk.params.bins)

#quick average of the profiles
mean_val = np.mean(val,axis=0)
std_val = np.std(val,axis=0)
np.savez('stack_mean_'+prof+'.npz',r=r[0],mean=mean_val,std=std_val)
#print "mean_val shape from profiles.py", np.shape(mean_val)
