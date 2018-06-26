import mpi4py.rc
import datetime
from pdict import *

def getparameters(filename):

    import params
    import globals

    dict=pdict()
    dict.read_from_file(filename)

    if 'basepath'      in dict: params.basepath      = dict['basepath']
    if 'serial'        in dict: params.serial        = dict['serial']
    if 'search_radius' in dict: params.search_radius = dict['search_radius']
    if 'lims'          in dict: params.lims          = dict['lims']
    if 'bins'          in dict: params.bins          = dict['bins']

    return params

def initialize(parameterfile):

    params = getparameters(parameterfile)

    if params.serial: mpi4py.rc.initialize = False
    from mpi4py import MPI

    if MPI.Is_initialized():
        params.comm     = MPI.COMM_WORLD
        params.rank     = params.comm.Get_rank()
        params.size     = params.comm.Get_size()
        params.parallel = True
    else:
        params.rank     = 0
        params.size     = 1
        params.parallel = False

    fmt       = '%H:%M:%S on %m/%d/%Y'
    timestamp = datetime.datetime.now().strftime(fmt)

    if(params.rank==0):
        print ''
        bar = 72*'-'
        print bar
        print 'Running on',params.size,'processor(s)'
        print 'Time:      '+timestamp
        print 'Directory: '+os.getcwd()
        print bar
        print ''

