#!/usr/bin/env python
import sys

import numpy             as np
import matplotlib.pyplot as plt
import illstack as istk

import mpi4py.rc
from decimal import Decimal

istk.init.initialize(sys.argv[1])
prof = str(sys.argv[2])  #This is getting specific profile (gas density, gas thermal pressure, dm)
mcenter=float(sys.argv[3])
ntile = 3 # controls tiling -- optimal value not yet clear

mlow=0.8*mcenter
mhigh=1.2*mcenter
mhmin = mlow / 1e10 # minimum mass in 1e10 Msun/h
mhmax = mhigh / 1e10 # maximum mass in 1e10 Msun/h
volweight = True # here we want density so value to bin is mass with a volume weight = true
scaled_radius=True
snap_num= int(sys.argv[4])
sim=str(sys.argv[5])

omegam=0.31
omegab=0.0486
omegadm =omegam-omegab
#Xh=0.76
gamma = 5./3.

mcenter=Decimal(mcenter)
mcenter='{:.2e}'.format(mcenter)
mcenter_power=np.log10(float(mcenter))
mcenter_power=str(int(mcenter_power))

if prof=='gasdens':
    print "Completing profiles for gasdens"
    part_type='gas'
    field_list = ['Coordinates','Masses']
    gas_particles = istk.io.getparticles(snap_num,part_type,field_list) #Change redshift
    posp = gas_particles['Coordinates'] #position, base unit ckpc/h 
    vals = gas_particles['Masses']   #units 1e10 Msol/h
elif prof=='dmdens':
    print "Completing profiles for dmdens"
    part_type='dm'
    # HARD CODED BOX SIZE 7.5e4 kpc/h
    part_massf=2.775e2*omegadm*(7.5e4)**3/1e10 # particle mass in 1e10 Msun/h
    field_list = ['Coordinates'] #base unit ckpc/h
    posp = istk.io.getparticles(snap_num,part_type,field_list) #Change redshift
    vals = posp[:,0]*0 + part_massf / np.shape(posp)[0]
    print 'setting dm particle mass to = ',vals[0]*1e10,'Msun/h'
elif prof=='gaspth':
    print "Completing profiles for gaspth"
    part_type='gas'
    field_list = ['Coordinates','Masses','InternalEnergy']
    gas_particles = istk.io.getparticles(snap_num,part_type,field_list) #Change redshift
    posp = gas_particles['Coordinates'] #base unit ckpc/h 
    vals = gas_particles['Masses']*gas_particles['InternalEnergy']*(gamma-1.) #unit 1e10Msol/h,(km/s)**2
else:
    print "Please enter an appropriate option for the profile"
    print "gasdens,dmdens,gaspth"


field_list = ['Group_M_Crit200','GroupPos','Group_R_Crit200','GroupFirstSub', 'GroupSFR','GroupMassType']
halos = istk.io.gethalos(snap_num,field_list)

GroupFirstSub=halos['GroupFirstSub']
posh = halos['GroupPos'] #kpc/h
mh   = halos['Group_M_Crit200']
rh   = halos['Group_R_Crit200'] #r200c, units ckpc/h
sfr  = halos['GroupSFR'] #Msol/yr
halo_mass= halos['GroupMassType']
mstar= halo_mass[:,4] #stellar mass, 1e10 Msol/h

print "orignal posh", np.shape(posh)
print "shape of groupfirstsub", np.shape(GroupFirstSub)

r, val, n, mh, rh, nprofs,GroupFirstSub,sfr,mstar = istk.cyprof.stackonhalos(posp,vals,posh,mh,rh,GroupFirstSub,sfr,mstar,ntile,volweight,mhmin, mhmax,scaled_radius)
print "nprofs", nprofs
r  =np.reshape(r,  (nprofs,istk.params.bins))
val=np.reshape(val,(nprofs,istk.params.bins)) 
n  =np.reshape(n,  (nprofs,istk.params.bins))

print 'shapes: ','r', np.shape(r),'val', np.shape(val),'n', np.shape(n),'mh', np.shape(mh)

np.savez(prof+'_scaled_'+sim+'_'+str(sys.argv[4])+'_'+mcenter_power+'.npz',r=r[0],val=val,n=n,mh=mh,rh=rh,nprofs=nprofs,nbins=istk.params.bins,GroupFirstSub=GroupFirstSub,sfr=sfr,mstar=mstar)
