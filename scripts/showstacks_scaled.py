import sys
import matplotlib.pyplot as plt
import numpy             as np
from scipy.interpolate import *

stacks_dens=np.load(sys.argv[1])
stacks_pres=np.load(sys.argv[2])
mlow   = float(sys.argv[3])
mhigh = float(sys.argv[4])
mhmin = mlow / 1e10
mhmax = mhigh / 1e10
type=str(sys.argv[5])

val_dens = stacks_dens['val'] #units h^2*1e10Msol/kpc^3
r   = stacks_pres['r'] #all the r for pres and dens are the same, and same for each profile
val_pres = stacks_pres['val']  #units (1e10Msol/h)(km/s)^2/volume 
bins=stacks_pres['nbins']
nprofs = stacks_dens['nprofs']
mh     = stacks_dens['mh'] #Halo mass  #now M_crit200, units 1e10 Msol/h
rh     = stacks_dens['rh'] #This should be r200c, units c*kpc/h

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

z=0.2
rhocrit_z=rhocrit*(omegam*(1+z)**3+omegalam)
G=6.67e-11*1.989e30/((3.086e19)**3) #G in units kpc^3/(Msol*s^2)

x_values=[]
normalized_pres=[] #size (27,20)
for i in np.arange(nprofs):
    #r200c=(3./4./np.pi/rhombar*mh[i]/200)**(1./3.)
    r200c=rh[i]
    #x_values.append(r/r200c)
    x_values.append(r) #our r values should now already be r/r200c
    P200c=200.*G*mh[i]*rhocrit_z*omegab/(omegam*2.*r200c)
    pressure=val_pres[i,:]
    pressure_divnorm=pressure/P200c
    normalized_pres.append(pressure_divnorm)
    
    if type == 'pres':
        plt.loglog(r,(val_pres[i,:]/P200c),c='b',alpha=0.3)
    elif type == 'dens':
        plt.loglog(r,val_dens[i,:]/rhocrit_z,c='r',alpha=0.2) #each are 20 long
    elif type == 'pres_av':
        pass
    elif type == 'dens_av':
        pass

mean_dens=np.median(val_dens, axis=0)
mean_xvals=np.mean(x_values, axis=0)

y_pres=np.median(normalized_pres,axis=0)
y_dens=mean_dens/rhocrit_z

#get errors
norm_dens=val_dens/rhocrit_z
errup_dens=np.percentile(norm_dens,84,axis=0)
errlow_dens=np.percentile(norm_dens,16,axis=0)
errors_dens=np.array([[errlow_dens,errup_dens]])
errors_dens=np.reshape(errors_dens,(2,bins))

norm_pres=normalized_pres
errup_pres=np.percentile(norm_pres,84,axis=0)
errlow_pres=np.percentile(norm_pres,16,axis=0)
errors_pres=np.array([[errlow_pres,errup_pres]])
errors_pres=np.reshape(errors_pres,(2,bins))

np.savetxt('/Users/emilymoser/Desktop/IllustrisTNG_profiles/Scaled/Density_12.txt', np.c_[mean_xvals, y_dens,errlow_dens,errup_dens])
np.savetxt('/Users/emilymoser/Desktop/IllustrisTNG_profiles/Scaled/Pressure_12.txt', np.c_[mean_xvals, y_pres,errlow_pres,errup_pres])

if type == 'pres' or type == 'pres_av':
    plt.loglog(mean_xvals, y_pres, color='black')
    plt.errorbar(mean_xvals, y_pres,yerr=errors_pres,zorder=3,fmt='none',color='black') 
    plt.ylabel(r'$P/P_{c}$',fontsize=15)
elif type == 'dens' or type == 'dens_av':
    plt.loglog(mean_xvals, y_dens, color='black')
    plt.errorbar(mean_xvals, y_dens,yerr=errors_dens,zorder=3,fmt='none',color='black')
    plt.ylabel(r'$\rho/\rho_{c}$',fontsize=15)
plt.xlabel(r'$r/r_{200c}$',fontsize=15) 
plt.title(r'z=0.2 halos with $ %s < M < %s M_\odot$' % (str(mlow), str(mhigh)))
plt.savefig('/Users/emilymoser/Desktop/IllustrisTNG_profiles/Scaled/Scaled_'+type+'_12.pdf', bbox_inches='tight')
plt.savefig('/Users/emilymoser/Desktop/IllustrisTNG_profiles/Scaled/Scaled_'+type+'_12.png', bbox_inches='tight')
plt.show()

plt.close()
