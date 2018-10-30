import sys

import matplotlib.pyplot as plt
import numpy             as np

from scipy.interpolate import *
from scipy.stats import binned_statistic
stacks_dens=np.load(sys.argv[1])
stacks_pres=np.load(sys.argv[2])

val_dens = stacks_dens['val']
r_pres   = stacks_pres['r']
val_pres = stacks_pres['val']


r      = stacks_dens['r']
nprofs = stacks_dens['nprofs']
mh     = stacks_dens['mh'] #Halo mass
rh     = stacks_dens['rh']


omegam = 0.3
omegab = 0.042
rhombar = 2.775e2*omegam # <rho_matter> in h^2 Msun/kpc^3
rhobbar = 2.775e2*omegab # <rho_baryon> in h^2 Msun/kpc^3      
rhodbar = rhombar - rhobbar

mh   *= 1e10 # convert halo mass to Msun 
val_dens *= 1e10 # convert to h^2 Msun/kpc^3
val_pres *= 1e10 # convert to h^2 Msun/kpc^3

for i in np.arange(nprofs):  #For every halo
    r200m = (3./4./np.pi/rhombar*mh[i]/200.)**(1./3.)

    plt.loglog(r/r200m,val_dens[i,:],c='r',alpha=0.3)
    #plt.loglog(r/r200m,val_pres[i,:],c='b',alpha=0.3)
#print "number of profiles", i+1

##This is the different part
#average the halo masses then use to calculate average r200m that
#all the profiles will use, instead of each having their own
#mean_pres=np.mean(val_pres, axis=0)
mean_dens=np.mean(val_dens, axis=0)
mean_mh=np.mean(mh)
r200m = (3./4./np.pi/rhombar*mean_mh/200.)**(1./3.)
x_values=r/r200m
plt.loglog(x_values, mean_dens, color='black', label='Average')
plt.xlabel(r'$r/r_{200m}$',fontsize=15)
plt.ylabel(r'Density',fontsize=15)
plt.legend()
plt.title(r'all halos with $M>10^{13} M_\odot/h$')
#plt.savefig('/Users/emilymoser/Desktop/Density_all.png')
#plt.savefig('/Users/emilymoser/Desktop/Density_all.pdf')
plt.show()

plt.close()
