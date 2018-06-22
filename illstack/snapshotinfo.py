import illustris_python as il
import numpy as np

def meanparticleproperty(snapshot_number,partType,field):
    basePath='/project/projectdirs/m3058/malvarez/illustris/ill_data'
    field_list = [field]
    gas_mass = il.snapshot.loadSubset(basePath,135,partType,fields=field_list)
    mean_property = np.mean(gas_mass,dtype='double')*1e10/0.704

    return mean_property







