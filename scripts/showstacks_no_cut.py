import sys
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy             as np
from scipy.interpolate import *
from decimal import Decimal
import illustris_python as il

print "running showstacks_no_cut.py"

stacks_dens=np.load(sys.argv[1])
stacks_pres=np.load(sys.argv[2])

mcenter=float(sys.argv[3])
mlow=0.8*mcenter
mhigh=1.2*mcenter


mhmin = mlow / 1e10
mhmax = mhigh / 1e10
snap=str(sys.argv[4])
sim=str(sys.argv[5])
type=str(sys.argv[6])

if sim =='ill':
    basepath='/global/cscratch1/sd/emm376/Illustris-3/output'
else:
    basepath='/global/cscratch1/sd/emm376/TNG100-3/output'

#mcenter=1.25*mlow
mcenter=Decimal(mcenter)
mcenter='{:.2e}'.format(mcenter)
mcenter_power=np.log10(float(mcenter))
mcenter_power=str(int(mcenter_power))

val_dens = stacks_dens['val'] #units h^2*1e10Msol/kpc^3
r   = stacks_pres['r'] #all the r for pres and dens are the same, and same for each profile
val_pres = stacks_pres['val']  #units (1e10Msol/h)(km/s)^2/volume 
bins=stacks_pres['nbins']
nprofs = stacks_dens['nprofs']
mh     = stacks_dens['mh'] #Halo mass  #now M_crit200, units 1e10 Msol/h
rh     = stacks_dens['rh'] #This should be r200c, units c*kpc/h
GroupFirstSub=stacks_dens['GroupFirstSub']
#print "max rh from showstacks", np.max(rh)

omegam = 0.31
omegalam=0.69
omegab = 0.0486
h=0.677
rhombar = 2.775e2*omegam # <rho_matter> in h^2 Msun/kpc^3
rhobbar = 2.775e2*omegab # <rho_baryon> in h^2 Msun/kpc^3      
rhodbar = rhombar - rhobbar
rhocrit=2.775e2

mh *= 1e10 # convert halo mass to Msun
val_dens *= 1e10 # convert to h^2 Msun/kpc^3 
val_pres *= 1e10 
val_pres /= (3.086e16*3.086e16)

#redshift
if snap == '050' or snap == '085':
    z=1.0
elif snap == '084' or snap =='120':
    z=0.2
elif snap == '064' or snap =='100':
    z=0.58

rhocrit_z=rhocrit*(omegam*(1+z)**3+omegalam)
G=6.67e-11*1.989e30/((3.086e19)**3) #G in units kpc^3/(Msol*s^2)

x_values=[]
unnormalized_pressure=[]
normalized_pres=[] #size (27,20)
for i in np.arange(nprofs):
    #r200c=(3./4./np.pi/rhombar*mh[i]/200)**(1./3.)
    r200c=rh[i]
    x_values.append(r) #our r values should now already be r/r200c
    P200c=200.*G*mh[i]*rhocrit_z*omegab/(omegam*2.*r200c)
    pressure=val_pres[i,:]
    unnormalized_pressure.append(pressure)
    pressure_divnorm=pressure/P200c
    normalized_pres.append(pressure_divnorm)
    
    if type == 'pres':
        plt.loglog(r,(val_pres[i,:]/P200c),c='darkgrey',alpha=0.3)
    elif type == 'dens':
        plt.loglog(r,val_dens[i,:]/rhocrit_z,c='darkgrey',alpha=0.2) #each are 20 long
    elif type == 'pres_av':
        pass
    elif type == 'dens_av':
        pass


mean_dens=np.median(val_dens, axis=0)
mean_xvals=np.mean(x_values, axis=0)

y_pres=np.median(normalized_pres,axis=0)
y_dens=mean_dens/rhocrit_z

#get unnormalized values
mean_unnorm_pres=np.median(unnormalized_pressure,axis=0)
#p200c=mean_unnorm_pres/y_pres
"""
#get errors
norm_dens=val_dens/rhocrit_z
errup_dens=np.percentile(norm_dens,84,axis=0)
errlow_dens=np.percentile(norm_dens,16,axis=0)
errors_dens=np.array([[errlow_dens,errup_dens]])
errors_dens=np.reshape(errors_dens,(2,bins))
"""
norm_pres=normalized_pres
errup_pres=np.percentile(norm_pres,84,axis=0)
errlow_pres=np.percentile(norm_pres,16,axis=0)
errors_pres=np.array([[errlow_pres,errup_pres]])
errors_pres=np.reshape(errors_pres,(2,bins))
"""
header_dens='x(r/r200c),rho/rhoc,rho(Msol/kpc^3),lower error(rho/rhoc),upper error(rho/rhoc)'
header_pres='x(r/r200c),P/P200c,P(Msol/kpc/s^2),lower error (P/P200c),upper error (P/P200c)' 
np.savetxt('/global/homes/e/emm376/Repositories/illstack/scripts/Results/Density_'+sim+'_'+snap+'_'+mcenter_power+'.txt', np.c_[mean_xvals,y_dens,mean_dens,errlow_dens,errup_dens],header=header_dens)
np.savetxt('/global/homes/e/emm376/Repositories/illstack/scripts/Results/Pressure_'+sim+'_'+snap+'_'+mcenter_power+'.txt', np.c_[mean_xvals,y_pres,mean_unnorm_pres,errlow_pres,errup_pres],header=header_pres)

if type == 'pres' or type == 'pres_av':
    plt.loglog(mean_xvals, y_pres, color='black', label='nprofs %s'%nprofs)
    #plt.errorbar(mean_xvals, y_pres,yerr=errors_pres,zorder=3,fmt='none',color='black')
    plt.fill_between(mean_xvals,y_pres+errup_pres, y_pres-errlow_pres,facecolor='grey',alpha=0.4)
    plt.ylabel(r'$P/P_{c}$',fontsize=15)
elif type == 'dens' or type == 'dens_av':
    plt.loglog(mean_xvals, y_dens, color='black',label='nprofs %s'%nprofs)
    #plt.errorbar(mean_xvals, y_dens,yerr=errors_dens,zorder=3,fmt='none',color='black')
    plt.fill_between(mean_xvals,y_dens+errup_dens, y_dens-errlow_dens,facecolor='grey',alpha=0.4)
    plt.ylabel(r'$\rho/\rho_{c}$',fontsize=15)
plt.legend()
plt.xlabel(r'$r/r_{200c}$',fontsize=15) 
plt.title(r'z= %s halos with $ %s < M < %s M_\odot$' % (str(z),str(mlow), str(mhigh)))
plt.savefig('/global/homes/e/emm376/Repositories/illstack/scripts/Results/July_10_'+type+'_'+sim+'_'+snap+'_'+mcenter_power+'.pdf', bbox_inches='tight')
plt.close()
"""
#merger trees
for c in np.arange(len(GroupFirstSub)):
    tree=il.sublink.loadTree(basepath,float(snap),GroupFirstSub[c],fields=['SubhaloStellarPhotometrics','SnapNum'],onlyMPB=True)
    mags=tree['SubhaloStellarPhotometrics']
    mag_g=mags[:,4]
    mag_r=mags[:,5]
    color=mag_g-mag_r
    plt.plot(tree['SnapNum'],color,'-')
plt.xlabel('Snapshot Number')
plt.ylabel('Color g-r')
plt.title('Color merger,nprofs = %s'%str(nprofs))
plt.savefig('/global/homes/e/emm376/Repositories/illstack/scripts/Results/July_11_'+sim+'_'+snap+'_merger_colors_total_'+mcenter_power+'.pdf', bbox_inches='tight')
plt.close()

for s in np.arange(len(GroupFirstSub)):
    tree=il.sublink.loadTree(basepath,float(snap),GroupFirstSub[s],fields=['SubhaloSFR','SnapNum'],onlyMPB=True)
    sfr=tree['SubhaloSFR']
    plt.plot(tree['SnapNum'],sfr,'-')
plt.xlabel('Snapshot Number')
plt.ylabel('Subhalo SFR (Msol/yr)')
plt.title('SFR merger,nprofs = %s'%str(nprofs))
plt.savefig('/global/homes/e/emm376/Repositories/illstack/scripts/Results/July_11_'+sim+'_'+snap+'_merger_sfr_total_'+mcenter_power+'.pdf', bbox_inches='tight')
plt.close()

for ms in np.arange(len(GroupFirstSub)):
    tree=il.sublink.loadTree(basepath,float(snap),GroupFirstSub[ms],fields=['SubhaloMassType','SnapNum'],onlyMPB=True)
    m=tree['SubhaloMassType']
    mstar=m[:,4]
    plt.plot(tree['SnapNum'],mstar,'-')
plt.xlabel('Snapshot Number')
plt.yscale('log')
plt.ylabel('Stellar Mass (1e10 Msol/h)')
plt.title('Mstar merger,nprofs = %s'%str(nprofs))
plt.savefig('/global/homes/e/emm376/Repositories/illstack/scripts/Results/July_11_'+sim+'_'+snap+'_merger_mstar_total_'+mcenter_power+'.pdf', bbox_inches='tight')
plt.close()
