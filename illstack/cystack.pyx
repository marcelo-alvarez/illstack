import healpy as hp
import numpy as np
cimport numpy as np

def some_fast_cython_stacking_routine(np.ndarray x, n):
    
    '''
    THIS IS A TEMPLATE FUNCTION FOR CYTHON IN ILLSTACK
    Parameters
    x[:]: an input array
    n: 	  some input integer

    Returns
    the_output[0:n-1]: an output array
    '''
    the_output = np.zeros(n)
    return the_output
