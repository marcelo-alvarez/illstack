import illustris_python as il
import numpy as np

basePath='/project/projectdirs/m3058/malvarez/illustris/ill_data'
fields = ['Masses']
gas_mass = il.snapshot.loadSubset(basePath,135,'gas',fields=fields)
print np.log10( np.mean(gas_mass,dtype='double')*1e10/0.704 )







