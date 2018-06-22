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
    
    nhalos=np.shape(xh)[0]
    
    ninhalos=0
    nphalo = np.zeros(nhalos)
    nhalos = 10
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
	
def stackonhalos(
    np.ndarray xp,np.ndarray yp,np.ndarray zp,
    np.ndarray xh,np.ndarray yh,np.ndarray zh,
    np.ndarray rh,np.ndarray mh, it, jt, kt,ntile):
      for i in np.arange(it):
        for j in np.arange(jt):
          for k in np.arange(kt):
             xpinc,ypinc,zpinc = stackonhalostile(xp,yp,zp,xh,yh,zh,rh,mh,
	                                          it,jt,kt)
             
		     

