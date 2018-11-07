import illustris_python as il
import numpy as np
import params

def getparticles(snapshot_number,partType,field_list):

    basePath=params.basepath
    particles = il.snapshot.loadSubset(basePath,snapshot_number,partType,fields=field_list)

    return particles

def gethalos(snapshot_number,field_list):
    
    basePath=params.basepath
    return il.groupcat.loadHalos(basePath,snapshot_number,fields=field_list)





