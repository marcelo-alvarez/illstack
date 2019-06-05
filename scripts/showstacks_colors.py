import sys
import matplotlib.pyplot as plt
import numpy             as np
from scipy.interpolate import *

#Right now this makes profiles separately for dens/pres, so go into loop and plotting to define
stacks_dens_red=np.load(sys.argv[1]) #these should be in base units 
stacks_pres_red=np.load(sys.argv[2])
stacks_dens_blue=np.load(sys.argv[3])
stacks_pres_blue=np.load(sys.argv[4])

mlow   = float(sys.argv[5])
mhigh = float(sys.argv[6])
mhmin = mlow / 1e10  # minimum mass in 1e10 Msun/h
mhmax = mhigh / 1e10 # maximum mass in 1e10 Msun/h

type=str(sys.argv[7])

val_dens_red = stacks_dens_red['val'] #units h^2*1e10Msol/kpc^3
val_dens_blue = stacks_dens_blue['val']
bins=stacks_pres_blue['nbins']

r_red   = stacks_pres_red['r'] #all the r for pres and dens are the same, and same for each profile
r_blue   = stacks_pres_blue['r']
val_pres_red = stacks_pres_red['val']  #units (1e10Msol/h)(km/s)^2/volume 
val_pres_blue = stacks_pres_blue['val']

nprofs_red = stacks_dens_red['nprofs']
mh_red     = stacks_dens_red['mh'] #Halo mass  #now M_crit200, units 1e10 Msol/h
rh_red     = stacks_dens_red['rh'] #This should be r200c, units c*kpc/h


nprofs_blue = stacks_dens_blue['nprofs']
mh_blue     = stacks_dens_blue['mh']
rh_blue     = stacks_dens_blue['rh']

#These are the IllustrisTNG values
omegam = 0.31
omegalam=0.69
omegab = 0.0486
h=0.677
rhombar = 2.775e2*omegam # <rho_matter> in h^2 Msun/kpc^3
rhobbar = 2.775e2*omegab # <rho_baryon> in h^2 Msun/kpc^3      
rhodbar = rhombar - rhobbar
rhocrit=2.775e2 #in units of h^2

mh_red *= 1e10 # convert halo mass to Msun
val_dens_red *= 1e10 # convert to h^2 Msun/kpc^3 
val_pres_red *= 1e10 # convert to h^2 Msun/kpc^3 #This should actually be Msol*km^2/s^2/kpc^3
#Convert val_pres into Msol/(s^2*kpc) (basically cancel the km**2 with kpc**2)
val_pres_red /= (3.086e16*3.086e16)

mh_blue *= 1e10
val_dens_blue *= 1e10
val_pres_blue *= 1e10
val_pres_blue /= (3.086e16*3.086e16)

#Normalize the y-axis by critical values                                        
z=0.2
rhocrit_z=rhocrit*(omegam*(1+z)**3+omegalam)
G=6.67e-11*1.989e30/((3.086e19)**3) #G in units kpc^3/(Msol*s^2)

x_values_red=[]
normalized_pres_red=[] #size (27,20)
for i in np.arange(nprofs_red):  #For every halo
    r200c=rh_red[i]
    x_values_red.append(r_red/r200c)
    P200c=200.*G*mh_red[i]*rhocrit_z*omegab/(omegam*2.*r200c)
    pressure=val_pres_red[i,:] #size 20
    pressure_divnorm=pressure/P200c  #size 20
    normalized_pres_red.append(pressure_divnorm)
    
    if type == 'pres':
        plt.loglog(r_red/r200c,(val_pres_red[i,:]/P200c),c='lightcoral',alpha=0.3)
    elif type == 'pres_av':
        pass
    elif type == 'dens':
        plt.loglog(r_red/r200c,val_dens_red[i,:]/rhocrit_z,c='lightcoral',alpha=0.3)
    elif type =='dens_av':
        pass

mean_dens_red=np.median(val_dens_red, axis=0)
mean_xvals_red=np.mean(x_values_red, axis=0)

#red errors
norm_dens_red=val_dens_red/rhocrit_z
errup_dens_red=np.percentile(norm_dens_red,84,axis=0)
errlow_dens_red=np.percentile(norm_dens_red,16,axis=0)
errors_dens_red=np.array([[errlow_dens_red,errup_dens_red]])
errors_dens_red=np.reshape(errors_dens_red,(2,bins))

norm_pres_red=normalized_pres_red
errup_pres_red=np.percentile(norm_pres_red,84,axis=0)
errlow_pres_red=np.percentile(norm_pres_red,16,axis=0)
errors_pres_red=np.array([[errlow_pres_red,errup_pres_red]])
errors_pres_red=np.reshape(errors_pres_red,(2,bins))

y_pres_red=np.median(normalized_pres_red,axis=0)
y_dens_red=mean_dens_red/rhocrit_z

#np.savetxt('/Users/emilymoser/Desktop/IllustrisTNG_profiles/Color_test/Mass_12/Density_red_12.txt', np.c_[mean_xvals_red, y_dens_red,errlow_dens_red,errup_dens_red])
#np.savetxt('/Users/emilymoser/Desktop/IllustrisTNG_profiles/Color_test/Mass_12/Pressure_red_12.txt', np.c_[mean_xvals_red, y_pres_red,errlow_pres_red,errup_pres_red])


#Now repeat everything for blue
x_values_blue=[]
normalized_pres_blue=[]
for i in np.arange(nprofs_blue):
    #r200c=(3./4./np.pi/rhombar*mh_blue[i]/200)**(1./3.)
    r200c=rh_blue[i]
    x_values_blue.append(r_blue/r200c)
    P200c=200.*G*mh_blue[i]*rhocrit_z*omegab/(omegam*2.*r200c)
    pressure=val_pres_blue[i,:]
    pressure_divnorm=pressure/P200c
    normalized_pres_blue.append(pressure_divnorm)
    
    if type == 'pres':
        plt.loglog(r_blue/r200c,(val_pres_blue[i,:]/P200c),c='cornflowerblue',alpha=0.3)
    elif type == 'pres_av':
        pass
    elif type == 'dens':
        plt.loglog(r_blue/r200c,val_dens_blue[i,:]/rhocrit_z,c='cornflowerblue',alpha=0.3)
    elif type == 'dens_av':
        pass

mean_dens_blue=np.median(val_dens_blue, axis=0)
mean_xvals_blue=np.mean(x_values_blue, axis=0)

#blue errors
norm_dens_blue=val_dens_blue/rhocrit_z
errup_dens_blue=np.percentile(norm_dens_blue,84,axis=0)
errlow_dens_blue=np.percentile(norm_dens_blue,16,axis=0)
errors_dens_blue=np.array([[errlow_dens_blue,errup_dens_blue]])
errors_dens_blue=np.reshape(errors_dens_blue,(2,bins))

norm_pres_blue=normalized_pres_blue
errup_pres_blue=np.percentile(norm_pres_blue,84,axis=0)
errlow_pres_blue=np.percentile(norm_pres_blue,16,axis=0)
errors_pres_blue=np.array([[errlow_pres_blue,errup_pres_red]])
errors_pres_blue=np.reshape(errors_pres_blue,(2,bins))

y_pres_blue=np.median(normalized_pres_blue,axis=0)
y_dens_blue=mean_dens_blue/rhocrit_z 

#np.savetxt('/Users/emilymoser/Desktop/IllustrisTNG_profiles/Color_test/Mass_12/Density_blue_12.txt', np.c_[mean_xvals_blue, y_dens_blue,errlow_dens_blue,errup_dens_blue])
#np.savetxt('/Users/emilymoser/Desktop/IllustrisTNG_profiles/Color_test/Mass_12/Pressure_blue_12.txt', np.c_[mean_xvals_blue, y_pres_blue,errlow_pres_blue,errup_pres_blue])

if type == 'pres' or type == 'pres_av':
    plt.loglog(mean_xvals_red, y_pres_red, color='red',label='red')
    plt.loglog(mean_xvals_blue, y_pres_blue, color='blue',label='blue')
    plt.errorbar(mean_xvals_red, y_pres_red,yerr=errors_pres_red,zorder=3,fmt='none',color='darkred')
    plt.errorbar(mean_xvals_blue, y_pres_blue,yerr=errors_pres_blue,zorder=3,fmt='none',color='darkblue')
    plt.ylabel(r'$P/P_{c}$',fontsize=15)

elif type == 'dens' or type == 'dens_av':
    plt.loglog(mean_xvals_red, y_dens_red, color='red', label='red')
    plt.loglog(mean_xvals_blue, y_dens_blue, color='blue', label='blue')
    plt.errorbar(mean_xvals_red, y_dens_red,yerr=errors_dens_red,zorder=3,fmt='none',color='darkred')
    plt.errorbar(mean_xvals_blue, y_dens_blue,yerr=errors_dens_blue,zorder=3,fmt='none',color='midnightblue')
    plt.ylabel(r'$\rho/\rho_{c}$',fontsize=15)

plt.xlabel(r'$r/r_{200c}$',fontsize=15)
plt.legend()
plt.title(r'z=0.2 halos with $ %s < M < %s M_\odot$' % (str(mlow), str(mhigh)))
#plt.savefig('/Users/emilymoser/Desktop/IllustrisTNG_profiles/Color_test/Mass_12/Density'+type+'_12.pdf', bbox_inches='tight')
#plt.savefig('/Users/emilymoser/Desktop/IllustrisTNG_profiles/Color_test/Mass_12/Density'+type+'_12_av.png', bbox_inches='tight')
plt.show()

plt.close()
