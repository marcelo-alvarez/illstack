import healpy as hp
import numpy as np
cimport numpy as np

search_radius = 2
ntile=8
def stackonhalostile(
    np.ndarray xp,np.ndarray yp,np.ndarray zp,
    np.ndarray xh,np.ndarray yh,np.ndarray zh,
    np.ndarray rh,np.ndarray mh, it, jt, kt,ntile):

    '''
    Parameters
	particles[nparticles][npartprops]
        halos[nphalos][nhaloprops]

    Returns
	profiles[:,nhalos]
    '''

    rpmax = rh.max()
    rbuff = rpmax * search_radius

    x1=xp.min(); x2=xp.max()
    y1=yp.min(); y2=yp.max()
    z1=zp.min(); z2=zp.max()

    dx=(x2-x1)/ntile; dy=(z2-z1)/ntile; dz=(z2-z1)/ntile;

    x1t=it*dx-rbuff; x2t=(it+1)*dx+rbuff
    y1t=jt*dy-rbuff; y2t=(jt+1)*dy+rbuff
    z1t=kt*dz-rbuff; z2t=(kt+1)*dz+rbuff

    dm =  [(xp>x1t) & (xp<x2t) & (yp>y1t) & (yp<y2t) & (zp>z1t) & (zp<z2t)]
    x1t+=rbuff;x2t-=rbuff;y1t+=rbuff;y2t-=rbuff;z1t+=rbuff;z2t-=rbuff;
    dmh = [(xh>x1t) & (xh<x2t) & (yh>y1t) & (yh<y2t) & (zh>z1t) & (zh<z2t)]

    xp=xp[dm];  yp=yp[dm];  zp=zp[dm]
    xh=xh[dmh]; yh=yh[dmh]; zh=zh[dmh]

    nhalos=np.shape(xh)[0]
    
    ninhalos=0
    nphalo = np.zeros(nhalos)
    xpin = []; ypin=[]; zpin=[]
    for ih in np.arange(nhalos):    	
        #        if ih%1 == 1000: print ih,'done out of',nhalos
        print ih,'done out of',nhalos
        rp = np.sqrt((xp-xh[ih])**2+(yp-yh[ih])**2+(zp-zh[ih])**2)
        dm = [np.abs(rp) < rh[ih]]
        
        xpin.extend(xp[dm])
        ypin.extend(yp[dm])
        zpin.extend(zp[dm])
        
        ninhalos += len(xp[dm])

    xpin = np.asarray(xpin)
    ypin = np.asarray(ypin)
    zpin = np.asarray(zpin)
    
    print 'ninhalos = ',ninhalos
    print 'size of xpin = ',len(xpin) 

    return xpin, ypin, zpin
	
def stackonhalos(
        np.ndarray xp,np.ndarray yp,np.ndarray zp,
        np.ndarray xh,np.ndarray yh,np.ndarray zh,
        np.ndarray rh,np.ndarray mh, ntile):

    pcen = np.empty((0),float)
    pval = np.empty((0),float)
    pnum = np.empty((0),float)
    
    for it in np.arange(ntile):
        for jt in np.arange(ntile):
            for kt in np.arange(ntile):
                print it, jt, kt
                pcenc, pvalc, pnumc = stackonhalostile(xp,yp,zp,xh,yh,zh,rh,mh,
                                                       it,jt,kt,ntile)
                pcen=np.append(pcen,pcenc)
                pval=np.append(pval,pvalc)
                pnum=np.append(pnum,pnumc)

    return pcen, pval, pnum
		     

