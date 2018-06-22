import illustris_python as il
import numpy as np
import params

def meanparticleproperty(snapshot_number,partType,field):
    basePath=params.basepath
    field_list = [field]
    gas_mass = il.snapshot.loadSubset(basePath,135,partType,fields=field_list)
    mean_property = np.mean(gas_mass,dtype='double')*1e10/0.704

    return mean_property







