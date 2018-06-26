import healpy as hp
import numpy as np
cimport numpy as np
import params

from illstack.CompHaloProperties import CompHaloProp
lims = [4,15.]
bins = 10
CHP = CompHaloProp(lims,bins)
search_radius = params.search_radius
box = 75000. # FIX THIS!!!!

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
        it, jt, kt,ntile,volweight):

    '''
    Parameters
	particles[nparticles][npartprops]
        halos[nphalos][nhaloprops]

    Returns
	profiles[:,nhalos]
    '''

    rpmax = rh.max()
    rbuff = rpmax * search_radius

    xp = posp[:,0]; yp = posp[:,1]; zp = posp[:,2]
    xh = posh[:,0]; yh = posh[:,1]; zh = posh[:,2]

    x1=xp.min(); x2=xp.max()
    y1=yp.min(); y2=yp.max()
    z1=zp.min(); z2=zp.max()

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

    if x1h>=x2h or y1h>=y2h or z1h>=z2h: 
        print 'error: too many tiles'
        print it,jt,kt,x1h,x2h,y1h,y2h,z1h,z2h
    
    dmp = [(xp>x1p) & (xp<x2p) & (yp>y1p) & (yp<y2p) & (zp>z1p) & (zp<z2p)]
    dmh = [(xh>x1h) & (xh<x2h) & (yh>y1h) & (yh<y2h) & (zh>z1h) & (zh<z2h)]

    xp=xp[dmp]; yp=yp[dmp]; zp=zp[dmp]; vals=vals[dmp];
    xh=xh[dmh]; yh=yh[dmh]; zh=zh[dmh]; mh=mh[dmh]; rh=rh[dmh]

    posp  = np.column_stack([xp,yp,zp])
    posh  = np.column_stack([xh,yh,zh])

    nhalos=np.shape(xh)[0]
    
    if params.rank==0:
        print 'after tiling nhalos = ',nhalos
        print it,jt,kt,'of',ntile,ntile,ntile,' done'
    
    ninhalos=0
    nphalo = np.zeros(nhalos)
    
    pcen = np.empty((0),float)
    pval = np.empty((0),float)
    pnum = np.empty((0),float)

    if nhalos == 0:
        return pcen, pval, pnum, mh

    weights = 1.0 + 0*xp

    for ih in np.arange(nhalos):    	

        pospc, valsc, weightsc = cull_and_center(posp,vals,weights,posh[ih],rh[ih])        

        pcenc, pvalc, pnumc = CHP.ComputeHaloProfile(pospc,valsc,weightsc,volweight=volweight)
        pcen = np.append(pcen,pcenc)
        pval = np.append(pval,pvalc)
        pnum = np.append(pnum,pnumc)
        if ih%100 == 0: print ih,'done out of',nhalos        

    return pcen,pval,pnum,mh
	
def stackonhalos(
        np.ndarray posp,
        np.ndarray vals,
        np.ndarray posh,
        np.ndarray   mh,
        np.ndarray   rh,
        ntile, volweight):
    
    pcen = np.empty((0),float)
    pval = np.empty((0),float)
    pnum = np.empty((0),float)
    mhpr = np.empty((0),float)

#    for it in np.arange(ntile):
#        for jt in np.arange(ntile):
#            for kt in np.arange(ntile):

    for it in np.arange(1):
        for jt in np.arange(1):
            for kt in np.arange(1):
                print it, jt, kt
                pcenc, pvalc, pnumc, mhc = stackonhalostile(posp,vals,posh,mh,rh,
                                                            it+2,jt+2,kt+2,ntile,volweight)
                pcen=np.append(pcen,pcenc)
                pval=np.append(pval,pvalc)
                pnum=np.append(pnum,pnumc)
                mhpr=np.append(mhpr,  mhc)

    return pcen, pval, pnum, mhpr
		     

