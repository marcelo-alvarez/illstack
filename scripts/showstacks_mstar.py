import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy             as np
from scipy.interpolate import *
from decimal import Decimal
import illustris_python as il

print "running showstacks_mstar.py"
cut_mstar=1e10
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

mcenter=Decimal(mcenter)
mcenter='{:.2e}'.format(mcenter)
mcenter_power=np.log10(float(mcenter))
mcenter_power=str(int(mcenter_power))

val_dens = stacks_dens['val']
bins=stacks_pres['nbins']

r   = stacks_pres['r']
val_pres = stacks_pres['val']

nprofs = stacks_dens['nprofs']
mh     = stacks_dens['mh']
rh     = stacks_dens['rh']
GroupFirstSub=stacks_dens['GroupFirstSub']
mstar=stacks_dens['mstar']

percentile_16=np.percentile(mstar,16)
percentile_30=np.percentile(mstar,30)
percentile_50=np.percentile(mstar,50)
percentile_70=np.percentile(mstar,70)
percentile_84=np.percentile(mstar,84)
n,bins_hist,patches=plt.hist(mstar,bins=15)
plt.vlines(percentile_16,0,len(mstar),color='red',label='16 percent %.2f'%percentile_16)
plt.vlines(percentile_30,0,len(mstar),color='purple',label='30 percent %.2f'%percentile_30)
plt.vlines(percentile_50,0,len(mstar),color='yellow',label='50 percent %.2f'%percentile_50)
plt.vlines(percentile_70,0,len(mstar),color='black',label='70 percent %.2f'%percentile_70)
plt.vlines(percentile_84,0,len(mstar),color='green',label='84 percent %.2f'%percentile_84)
plt.xlabel('Mstar 1e10 Msol/h')
plt.title('Halo Stellar Mass, nprofs = %s'%str(nprofs))
plt.legend()
plt.ylim(0,max(n)+2)
#plt.savefig('/global/homes/e/emm376/Repositories/illstack/scripts/Results/July_9_'+type+'_'+sim+'_'+snap+'_mstar_histogram_'+mcenter_power+'.pdf', bbox_inches='tight')
plt.close()
#These are the IllustrisTNG values
omegam = 0.31
omegalam=0.69
omegab = 0.0486
h=0.677
rhombar = 2.775e2*omegam # <rho_matter> in h^2 Msun/kpc^3
rhobbar = 2.775e2*omegab # <rho_baryon> in h^2 Msun/kpc^3      
rhodbar = rhombar - rhobbar
rhocrit=2.775e2 #in units of h^2

mh *= 1e10
val_dens *= 1e10 
val_pres *= 1e10
val_pres /= (3.086e16*3.086e16)
mstar *= 1e10
#get redshift from snapshot number
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
normalized_pres=[]
for i in np.arange(nprofs):
    r200c=rh[i]
    x_values.append(r) #if scaled_radius is True
    P200c=200.*G*mh[i]*rhocrit_z*omegab/(omegam*2.*r200c)
    pressure=val_pres[i,:] #size 20
    unnormalized_pressure.append(pressure)
    pressure_divnorm=pressure/P200c  #size 20
    normalized_pres.append(pressure_divnorm)
    """
    if type == 'pres':
        plt.loglog(r_red,(val_pres_red[i,:]/P200c),c='lightcoral',alpha=0.3) #if scaled radius
    elif type == 'pres_av':
        pass
    elif type == 'dens':
        plt.loglog(r_red,val_dens_red[i,:]/rhocrit_z,c='lightcoral',alpha=0.3)
    elif type =='dens_av':
        pass
    """

percentile_high=np.percentile(mstar,70)
percentile_low=np.percentile(mstar,30)

#if we want to cut 
#idx_red=np.where(mstar >= cut_mstar)
#idx_blue=np.where(mstar < cut_mstar)
#if we want to split by percentile
idx_red=np.where(mstar >= percentile_high)
idx_blue=np.where(mstar <= percentile_low)

idx_red=np.array(idx_red[0])
idx_blue=np.array(idx_blue[0])

#print "mstar from showstacks", mstar
#print "mh", mh
nprofs_red=len(idx_red)
nprofs_blue=len(idx_blue)

print "nprofs high", nprofs_red
print "nprofs low", nprofs_blue
x_values=np.array(x_values)
normalized_pres=np.array(normalized_pres)
unnormalized_pressure=np.array(unnormalized_pressure)
val_dens=np.array(val_dens)

#red is high
if nprofs_red == 0:
    print "from showstacks_mstar.py we have 0 profiles with high mstar"
else:
    GroupFirstSub_red=GroupFirstSub[idx_red]
    mstar_red=mstar[idx_red]
    #print "mstar red", mstar_red
    x_values_red=x_values[idx_red]
    normalized_pres_red=normalized_pres[idx_red]
    unnormalized_pressure_red=unnormalized_pressure[idx_red]
    val_dens_red=val_dens[idx_red]
    mean_dens_red=np.median(val_dens_red, axis=0)
    mean_xvals_red=np.mean(x_values_red, axis=0)
    mean_unnorm_pres_red=np.median(unnormalized_pressure_red,axis=0)

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

    header_dens='x(r/r200c),rho/rhoc,rho(Msol/kpc^3),lower error(rho/rhoc),upper error(rho/rhoc)'
    header_pres='x(r/r200c),P/P200c,P(Msol/kpc/s^2),lower error (P/P200c),upper error (P/P200c)'
    #np.savetxt('/global/homes/e/emm376/Repositories/illstack/scripts/Results/Density_'+sim+'_'+snap+'_mstar_high_'+mcenter_power+'.txt', np.c_[mean_xvals_red,y_dens_red,mean_dens_red,errlow_dens_red,errup_dens_red],header=header_dens)
    #np.savetxt('/global/homes/e/emm376/Repositories/illstack/scripts/Results/Pressure_'+sim+'_'+snap+'_mstar_high_'+mcenter_power+'.txt', np.c_[mean_xvals_red,y_pres_red,mean_unnorm_pres_red,errlow_pres_red,errup_pres_red],header=header_pres)

#blue is low
if nprofs_blue == 0:
    print "from showstacks, we have 0 profiles with low mstar"
else:
    GroupFirstSub_blue=GroupFirstSub[idx_blue]
    mstar_blue=mstar[idx_blue]
    #print "mstar blue", mstar_blue
    x_values_blue=x_values[idx_blue]
    normalized_pres_blue=normalized_pres[idx_blue]
    unnormalized_pressure_blue=unnormalized_pressure[idx_blue]
    val_dens_blue=val_dens[idx_blue]
    mean_dens_blue=np.median(val_dens_blue, axis=0)
    mean_xvals_blue=np.mean(x_values_blue, axis=0)
    mean_unnorm_pres_blue=np.median(unnormalized_pressure_blue,axis=0)
    norm_dens_blue=val_dens_blue/rhocrit_z
    errup_dens_blue=np.percentile(norm_dens_blue,84,axis=0)
    errlow_dens_blue=np.percentile(norm_dens_blue,16,axis=0)
    errors_dens_blue=np.array([[errlow_dens_blue,errup_dens_blue]])
    errors_dens_blue=np.reshape(errors_dens_blue,(2,bins))

    norm_pres_blue=normalized_pres_blue
    errup_pres_blue=np.percentile(norm_pres_blue,84,axis=0)
    errlow_pres_blue=np.percentile(norm_pres_blue,16,axis=0)
    errors_pres_blue=np.array([[errlow_pres_blue,errup_pres_blue]])
    errors_pres_blue=np.reshape(errors_pres_blue,(2,bins))
    
    y_pres_blue=np.median(normalized_pres_blue,axis=0)
    y_dens_blue=mean_dens_blue/rhocrit_z
    header_dens='x(r/r200c),rho/rhoc,rho(Msol/kpc^3),lower error(rho/rhoc),upper error(rho/rhoc)'
    header_pres='x(r/r200c),P/P200c,P(Msol/kpc/s^2),lower error (P/P200c),upper error (P/P200c)'
    #np.savetxt('/global/homes/e/emm376/Repositories/illstack/scripts/Results/Density_'+sim+'_'+snap+'_mstar_low_'+mcenter_power+'.txt', np.c_[mean_xvals_blue,y_dens_blue,mean_dens_blue,errlow_dens_blue,errup_dens_blue],header=header_dens)
    #np.savetxt('/global/homes/e/emm376/Repositories/illstack/scripts/Results/Pressure_'+sim+'_'+snap+'_mstar_low_'+mcenter_power+'.txt', np.c_[mean_xvals_blue,y_pres_blue,mean_unnorm_pres_blue,errlow_pres_blue,errup_pres_blue],header=header_pres)
"""
if type == 'pres' or type == 'pres_av':
    if nprofs_red==0:
        pass
    else:
        plt.loglog(mean_xvals_red, y_pres_red, color='red',label='high nprofs %s'%str(nprofs_red))
        plt.fill_between(mean_xvals_red,y_pres_red+errup_pres_red,y_pres_red-errlow_pres_red,facecolor='lightcoral',alpha=0.5)
        plt.ylabel(r'$P/P_{c}$',fontsize=15)
    if nprofs_blue ==0:
        pass
    else:
        plt.loglog(mean_xvals_blue, y_pres_blue, color='blue',label='low nprofs %s'%str(nprofs_blue))
        plt.fill_between(mean_xvals_blue,y_pres_blue+errup_pres_blue, y_pres_blue-errlow_pres_blue,facecolor='royalblue',alpha=0.4)

elif type == 'dens' or type == 'dens_av':
    if nprofs_red==0:
        pass
    else:
        plt.loglog(mean_xvals_red, y_dens_red, color='red', label='high nprofs %s'%str(nprofs_red))
        plt.fill_between(mean_xvals_red,y_dens_red+errup_dens_red,y_dens_red-errlow_dens_red,facecolor='lightcoral',alpha=0.5)
        plt.ylabel(r'$\rho/\rho_{c}$',fontsize=15)
    if nprofs_blue == 0:
        pass
    else:
        plt.loglog(mean_xvals_blue, y_dens_blue, color='blue', label='low nprofs %s'%str(nprofs_blue))
        plt.fill_between(mean_xvals_blue,y_dens_blue+errup_dens_blue, y_dens_blue-errlow_dens_blue,facecolor='royalblue',alpha=0.4)
plt.xlabel(r'$r/r_{200c}$',fontsize=15)
plt.legend()
#plt.title(r'z=%s halos with $ %s < M < %s M_\odot$, mstar cut %s' % (str(z),str(mlow), str(mhigh),str(cut_mstar)))
plt.title(r'z=%s halos with $ %s < M < %s M_\odot$, 30 percent' % (str(z),str(mlow), str(mhigh)))
plt.savefig('/global/homes/e/emm376/Repositories/illstack/scripts/Results/July_8_'+type+'_'+sim+'_'+snap+'_mstar_percentile_'+mcenter_power+'.pdf', bbox_inches='tight')

plt.close()
"""
if nprofs_red == 0:
    pass
else:
    for r in np.arange(len(GroupFirstSub_red)):
        tree=il.sublink.loadTree(basepath,float(snap),GroupFirstSub_red[r],fields=['SubhaloMassType','SnapNum'],onlyMPB=True)
        m=tree['SubhaloMassType']
        mstar=m[:,4]
        plt.plot(tree['SnapNum'],mstar,'-',color='r',alpha=0.5)
if nprofs_blue ==0:
    pass
else:
    for b in np.arange(len(GroupFirstSub_blue)):
        tree=il.sublink.loadTree(basepath,float(snap),GroupFirstSub_blue[b],fields=['SubhaloMassType','SnapNum'],onlyMPB=True)
        m=tree['SubhaloMassType']
        mstar=m[:,4]
        plt.plot(tree['SnapNum'],mstar,'-',color='b',alpha=0.4)

plt.xlabel('Snapshot Number')
plt.yscale('log')
plt.ylabel('Stellar Mass (1e10 Msol/h)')
red_patch=mpatches.Patch(color='red',label='nprofs red %s'%str(nprofs_red))
blue_patch=mpatches.Patch(color='blue', label='nprofs blue %s'%str(nprofs_blue))
plt.legend(handles=[red_patch,blue_patch])
plt.title('Mstar merger,nprofs = %s'%str(nprofs))
plt.savefig('/global/homes/e/emm376/Repositories/illstack/scripts/Results/July_11_'+sim+'_'+snap+'_merger_mstar_cut_'+mcenter_power+'.pdf', bbox_inches='tight')
plt.close()
