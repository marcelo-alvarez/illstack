import healpy as hp
import numpy as np
cimport numpy as np
import params

from illstack.CompHaloProperties import CompHaloProp
search_radius = params.search_radius
box = 75000. # NEED TO FIX THIS!!!!

def cull_and_center(np.ndarray posp, np.ndarray vals, np.ndarray weights, np.ndarray posh, rh):

    xp = posp[:,0]-posh[0]; yp=posp[:,1]-posh[1]; zp=posp[:,2]-posh[2]
    r = np.sqrt(
        xp**2+yp**2+zp**2
        ) 
    dm = [r < search_radius * rh]

    xp=xp[dm];yp=yp[dm];zp=zp[dm];vals=vals[dm];weights=weights[dm]
    posp=np.column_stack([xp,yp,zp])

    return posp,vals,weights

def stackonhalostile(
        np.ndarray posp,
        np.ndarray   vals,
        np.ndarray posh,
        np.ndarray   mh,
        np.ndarray   rh,
        it, jt, kt,ntile,volweight,mhmin):

    '''
    Parameters
	particles[nparticles][npartprops]
        halos[nphalos][nhaloprops]

    Returns
	profiles[:,nhalos]
    '''

    CHP = CompHaloProp(params.lims,params.bins)

    rpmax = rh.max()
    rbuff = rpmax * search_radius

    xp = posp[:,0]; yp = posp[:,1]; zp = posp[:,2]
    xh = posh[:,0]; yh = posh[:,1]; zh = posh[:,2]

    x1=0.; x2=box; y1=0.; y2=box; z1=0.; z2=box;
    dx=(x2-x1)/ntile; dy=(y2-y1)/ntile; dz=(z2-z1)/ntile;

    x1p=it*dx; x2p=(it+1)*dx
    y1p=jt*dy; y2p=(jt+1)*dy
    z1p=kt*dz; z2p=(kt+1)*dz

    # select halos sufficiently far from box edges so as not to have to do periodic wrapping
    x1h=max(x1p,rbuff); x2h=min(x2p,box-rbuff)
    y1h=max(y1p,rbuff); y2h=min(y2p,box-rbuff)
    z1h=max(z1p,rbuff); z2h=min(z2p,box-rbuff)

    x1p=x1h-rbuff; x2p=x2h+rbuff
    y1p=y1h-rbuff; y2p=y2h+rbuff
    z1p=z1h-rbuff; z2p=z2h+rbuff

    dmp = [(xp>x1p) & (xp<x2p) & (yp>y1p) & (yp<y2p) & (zp>z1p) & (zp<z2p)]
    dmh = [(xh>x1h) & (xh<x2h) & (yh>y1h) & (yh<y2h) & (zh>z1h) & (zh<z2h) & (mh>mhmin)]

    xp=xp[dmp]; yp=yp[dmp]; zp=zp[dmp]; vals=vals[dmp];
    xh=xh[dmh]; yh=yh[dmh]; zh=zh[dmh]; mh=mh[dmh]; rh=rh[dmh]

    posp  = np.column_stack([xp,yp,zp])
    posh  = np.column_stack([xh,yh,zh])

    pcen = np.empty((0),float)
    pval = np.empty((0),float)
    pnum = np.empty((0),float)

    nhalos=np.shape(xh)[0]
    if params.rank==0:
        print it*ntile**2+jt*ntile+kt+1,'of',ntile**3,'done, nhalos =',nhalos
    	
    if nhalos == 0:
        return pcen, pval, pnum, mh, nhalos
    
    ninhalos=0
    nphalo = np.zeros(nhalos)
    
    weights = 1.0 + 0*xp

    for ih in np.arange(nhalos):    	

        pospc, valsc, weightsc = cull_and_center(posp,vals,weights,posh[ih],rh[ih])        

        pcenc, pvalc, pnumc = CHP.ComputeHaloProfile(pospc,valsc,weightsc,volweight=volweight)
        pcen = np.append(pcen,pcenc)
        pval = np.append(pval,pvalc)
        pnum = np.append(pnum,pnumc)

    return pcen,pval,pnum,mh,nhalos
	
def stackonhalos(
        np.ndarray posp,
        np.ndarray vals,
        np.ndarray posh,
        np.ndarray   mh,
        np.ndarray   rh,
        ntile, volweight,mhmin):
    
    pcen = np.empty((0),float)
    pval = np.empty((0),float)
    pnum = np.empty((0),float)
    mhpr = np.empty((0),float)

    nhalos=0
    for it in np.arange(ntile):
        for jt in np.arange(ntile):
            for kt in np.arange(ntile):

                pcenc, pvalc, pnumc, mhc, nhalosc = stackonhalostile(posp,vals,posh,mh,rh,
                                                                     it,jt,kt,ntile,volweight,mhmin)
                pcen=np.append(pcen,pcenc)
                pval=np.append(pval,pvalc)
                pnum=np.append(pnum,pnumc)
                mhpr=np.append(mhpr,  mhc)

                nhalos += nhalosc
                
    return pcen, pval, pnum, mhpr,nhalos
		     

