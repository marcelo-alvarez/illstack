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
#scaled_radius=True
snap_num= int(sys.argv[5])

omegam=0.31
omegab=0.0486
omegadm =omegam-omegab
#Xh=0.76
gamma = 5./3.

if prof=='gasdens':
    part_type='gas'
    field_list = ['Coordinates','Masses']
    gas_particles = istk.io.getparticles(snap_num,part_type,field_list) #Change redshift
    posp = gas_particles['Coordinates'] #position, base unit ckpc/h 
    vals = gas_particles['Masses']   #units 1e10 Msol/h
elif prof=='dmdens':
    part_type='dm'
    # HARD CODED BOX SIZE 7.5e4 kpc/h  ##also change something here?
    part_massf=2.775e2*omegadm*(7.5e4)**3/1e10 # particle mass in 1e10 Msun/h
    field_list = ['Coordinates'] #base unit ckpc/h
    posp = istk.io.getparticles(snap_num,part_type,field_list) #Change redshift
    vals = posp[:,0]*0 + part_massf / np.shape(posp)[0]
    print 'setting dm particle mass to = ',vals[0]*1e10,'Msun/h'
elif prof=='gaspth':
    part_type='gas'
    field_list = ['Coordinates','Masses','InternalEnergy']
    gas_particles = istk.io.getparticles(snap_num,part_type,field_list) #Change redshift
    posp = gas_particles['Coordinates'] #base unit ckpc/h 
    vals = gas_particles['Masses']*gas_particles['InternalEnergy']*(gamma-1.) #unit 1e10Msol/h,(km/s)**2
else:
    print "Please enter an appropriate option for the profile"
    print "gasdens,dmdens,gaspth"


field_list_halos = ['Group_M_Crit200','GroupPos','Group_R_Crit200','GroupFirstSub']
field_list_subhalos = ['SubhaloStellarPhotometrics']
halos = istk.io.gethalos(snap_num,field_list_halos)
subhalos=istk.io.getsubhalos(snap_num, field_list_subhalos)

primary_subhalo_indices=halos['GroupFirstSub']
print "number of halos", np.shape(primary_subhalo_indices)

mag_r=subhalos[:,5]
mag_g=subhalos[:,4]
color_gr=mag_g-mag_r
print "original color_gr", np.shape(color_gr)
#get the colors for the most massive subhalos of each halo
color_gr=color_gr[primary_subhalo_indices]
print "color gr after getting primaru subhalo indices", np.shape(color_gr)

idx_red=np.where(color_gr >= 0.5)
print "idx red", np.shape(idx_red)
idx_blue=np.where(color_gr < 0.5)
print "idx blue", np.shape(idx_blue)


posh = halos['GroupPos']
print "posh", np.shape(posh)
mh   = halos['Group_M_Crit200'] #1e10 Msol/h
print "mh", np.shape(mh)
rh   = halos['Group_R_Crit200']
print "rh", np.shape(rh)

posh_red=posh[idx_red]
print "posh_red", np.shape(posh_red)
posh_blue=posh[idx_blue]
print "posh blue", np.shape(posh_blue)
mh_red=mh[idx_red]
mh_blue=mh[idx_blue]
rh_red=rh[idx_red]
rh_blue=rh[idx_blue]

r_red, val_red, n_red, mh_red, rh_red, nprofs_red = istk.cyprof.stackonhalos(posp,vals,posh_red,mh_red,rh_red,ntile,volweight,mhmin, mhmax)
print "nprofs red", nprofs_red
r_red  =np.reshape(r_red,  (nprofs_red,istk.params.bins))
val_red=np.reshape(val_red,(nprofs_red,istk.params.bins)) 
n_red  =np.reshape(n_red,  (nprofs_red,istk.params.bins))

r_blue, val_blue, n_blue, mh_blue, rh_blue, nprofs_blue = istk.cyprof.stackonhalos(posp,vals,posh_blue,mh_blue,rh_blue,ntile,volweight,mhmin, mhmax)
print "nprofs blue", nprofs_blue
r_blue  =np.reshape(r_blue,  (nprofs_blue,istk.params.bins))
val_blue=np.reshape(val_blue,(nprofs_blue,istk.params.bins))
n_blue  =np.reshape(n_blue,  (nprofs_blue,istk.params.bins))

print 'red shapes: ','r_red', np.shape(r_red),'val_red', np.shape(val_red),'n_red', np.shape(n_red),'mh_red', np.shape(mh_red)
print 'blue shapes: ','r_blue', np.shape(r_blue),'val_blue', np.shape(val_blue),'n_blue', np.shape(n_blue),'mh_blue', np.shape(mh_blue)

np.savez('stack_'+prof+'_red_11.npz',r=r_red[0],val=val_red,n=n_red,mh=mh_red,rh=rh_red,nprofs=nprofs_red,nbins=istk.params.bins)
np.savez('stack_'+prof+'_blue_11.npz',r=r_blue[0],val=val_blue,n=n_blue,mh=mh_blue,rh=rh_blue,nprofs=nprofs_blue,nbins=istk.params.bins)
