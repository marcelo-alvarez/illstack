#!/usr/bin/env python
import sys
import numpy             as np
import matplotlib.pyplot as plt
import illstack as istk
import mpi4py.rc
from decimal import Decimal
import params_tng

istk.init.initialize(sys.argv[1])

prof = str(sys.argv[2])
mcenter=float(sys.argv[3])
snap_num= int(sys.argv[4])
sim=str(sys.argv[5])

ntile = 3 # controls tiling -- optimal value not yet clear

mass_option = params_tng.mass_option
mass_kind=params_tng.mass_kind #this is either 'halo' or 'stellar' 
if mass_option == 1:
    mlow=params_tng.mass_low
    mhigh=params_tng.mass_high
    print "Proceeding with mass option 1 for %s mass: upper/lower mass range, %s < M < %s"%(mass_kind,str(mlow),str(mhigh))
elif mass_option == 2:
    percentage = params_tng.mass_option_percent
    mlow = (1.0-percentage)*mcenter
    mhigh = (1.0+percentage)*mcenter
    print "Proceeding with mass option 2 for %s mass: central mass and %s percent ranges,%s < M < %s"%(mass_kind,100*percentage,str(mlow),str(mhigh))

mhmin = mlow /1e10 # minimum mass in 1e10 Msun/h
mhmax = mhigh /1e10 # maximum mass in 1e10 Msun/h

volweight = params_tng.volweight
scaled_radius=params_tng.scaled_radius

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

r, val, n, mh, rh, nprofs,GroupFirstSub,sfr,mstar = istk.cyprof.stackonhalos(posp,vals,posh,mh,rh,GroupFirstSub,sfr,mstar,ntile,volweight,mhmin, mhmax,scaled_radius,mass_kind)
print "nprofs", nprofs
r  =np.reshape(r,  (nprofs,istk.params.bins))
val=np.reshape(val,(nprofs,istk.params.bins)) 
n  =np.reshape(n,  (nprofs,istk.params.bins))

print 'shapes: ','r', np.shape(r),'val', np.shape(val),'n', np.shape(n),'mh', np.shape(mh)

#np.savez(prof+'_scaled_'+sim+'_'+str(sys.argv[4])+'_'+mcenter_power+'.npz',r=r[0],val=val,n=n,mh=mh,rh=rh,nprofs=nprofs,nbins=istk.params.bins,GroupFirstSub=GroupFirstSub,sfr=sfr,mstar=mstar)
np.savez(prof+'_unscaled_'+sim+'_'+str(sys.argv[4])+'_mstar_total_test100.npz',r=r[0],val=val,n=n,mh=mh,rh=rh,nprofs=nprofs,nbins=istk.params.bins,GroupFirstSub=GroupFirstSub,sfr=sfr,mstar=mstar)

#np.savez('test.npz',r=r[0],val=val,n=n,mh=mh,rh=rh,nprofs=nprofs,nbins=istk.params.bins,GroupFirstSub=GroupFirstSub,sfr=sfr,mstar=mstar)
