import numpy as np
import subprocess

#sim=np.array(['tng','ill'])
sim=np.array(['tng'])
mass=np.array([1e12])
mcenter_power=np.array([12])
#cut=np.array(['no_cut','color','sfr','mstar'])
cut=np.array(['no_cut'])

#Change the scaled radius and periodic boundary conditions in profiles_mult.py
prof1='dmdens'
prof2='gasdens'
prof3='gaspth'
type1='pres'
type2='dens'
type3='pres_av'
type4='dens_av'

for k in np.arange(len(sim)):
    if sim[k]=='ill':
        #snap=np.array(['085',100,120])
        snap=np.array(['120'])
    else:
        #snap=np.array(['084','064','050'])
        snap=np.array(['084'])
    for j in np.arange(len(snap)):
        for i in np.arange(len(mass)):
            f=open('getprof_temp_profiles.sh','w')
            print >>f, '#!/bin/bash'
            #print >> f, 'python', '../scripts/profiles_mult.py', 'istk-params_scaled.txt', prof1, str(mass[i]), snap[j], sim[k]  
            #print >> f, 'python', '../scripts/profiles_mult.py', 'istk-params_'+sim[k]+'.txt', prof2, str(mass[i]), snap[j], sim[k]
            print >> f, 'python', '../scripts/profiles_mult.py', 'istk-params_'+sim[k]+'.txt', prof3, str(mass[i]), snap[j], sim[k]
            f.close()
            print "*******************************************************************************************"
            print "profile shell generated for sim = %s, mass= %s, snap= %s"%(sim[k],str(mcenter_power[i]),str(snap[j]))
            print "now running the shell script"
            print "*******************************************************************************************"
            subprocess.call(['./getprof_temp_profiles.sh'],shell=True)                                                                           
            for m in np.arange(len(cut)):
                f=open('getprof_temp_plotting.sh','w')
                print >>f, '#!/bin/bash'
                print >> f, 'python', '../scripts/showstacks_'+cut[m]+'.py', prof2+'_scaled_'+sim[k]+'_'+str(snap[j])+'_'+str(mcenter_power[i])+'.npz', prof3+'_scaled_'+sim[k]+'_'+str(snap[j])+'_'+str(mcenter_power[i])+'.npz', str(mass[i]), snap[j], sim[k],type1                                                                                                                                     
                print >> f, 'python', '../scripts/showstacks_'+cut[m]+'.py', prof2+'_scaled_'+sim[k]+'_'+str(snap[j])+'_'+str(mcenter_power[i])+'.npz', prof3+'_scaled_'+sim[k]+'_'+str(snap[j])+'_'+str(mcenter_power[i])+'.npz',str(mass[i]), snap[j], sim[k],type2  
                print >> f, 'python', '../scripts/showstacks_'+cut[m]+'.py', prof2+'_scaled_'+sim[k]+'_'+str(snap[j])+'_'+str(mcenter_power[i])+'.npz', prof3+'_scaled_'+sim[k]+'_'+str(snap[j])+'_'+str(mcenter_power[i])+'.npz',str(mass[i]), snap[j], sim[k],type3
                print >> f, 'python', '../scripts/showstacks_'+cut[m]+'.py', prof2+'_scaled_'+sim[k]+'_'+str(snap[j])+'_'+str(mcenter_power[i])+'.npz', prof3+'_scaled_'+sim[k]+'_'+str(snap[j])+'_'+str(mcenter_power[i])+'.npz',str(mass[i]), snap[j], sim[k],type4

                f.close()
                print "-------------------------------------------------------------------------------------------"
                print "plotting shell generated for sim = %s, mass= %s, snap= %s"%(sim[k],str(mcenter_power[i]),str(snap[j]))
                print "now running the shell script"
                print "-------------------------------------------------------------------------------------------"
                #subprocess.call(['./getprof_temp_plotting.sh'],shell=True)
