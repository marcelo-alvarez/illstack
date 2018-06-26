import sys

import matplotlib.pyplot as plt
import numpy             as np

from scipy.interpolate import *

stacksd=np.load(sys.argv[1])
stacksg=np.load(sys.argv[2])

vald = stacksd['val']
rg   = stacksg['r']
valg = stacksg['val']

r      = stacksd['r']
nprofs = stacksd['nprofs']
mh     = stacksd['mh']

omegam = 0.3
omegab = 0.042
rhombar = 2.775e2*omegam # <rho_matter> in h^2 Msun/kpc^3
rhobbar = 2.775e2*omegab # <rho_baryon> in h^2 Msun/kpc^3
rhodbar = rhombar - rhobbar

mh   *= 1e10 # convert halo mass to Msun
vald *= 1e10 # convert to h^2 Msun/kpc^3
valg *= 1e10 # convert to h^2 Msun/kpc^3

rhogmean=0.
rhodmean=0.
for i in np.arange(nprofs):

    r200m = (3./4./np.pi/rhombar*mh[i]/200.)**(1./3.)
    rhodcur = vald[i,:]/rhodbar
    rhogcur = valg[i,:]/rhobbar
    xcur    = r/r200m

    plt.loglog(r/r200m,rhodcur,c='r',alpha=0.3)
    plt.loglog(r/r200m,rhogcur,c='b',alpha=0.3)

    rhodf = interp1d(xcur,rhodcur)
    rhogf = interp1d(xcur,rhogcur)    

plt.gca().set_xlabel(r'$r/r_{200m}$',fontsize=20)
plt.gca().set_ylabel(r'$\rho/\langle\rho\rangle$',fontsize=20)
plt.gca().set_title(r'all halos with $M>10^{13} M_\odot/h$')

plt.show()
