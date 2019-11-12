import healpy as hp
import numpy as np
cimport numpy as np
import params

from illstack.CompHaloProperties import CompHaloProp
search_radius = params.search_radius
box = 75000. # NEED TO FIX THIS!!!!

def periodic_bcs(np.ndarray posp,np.ndarray posh):
    
    xp = posp[:,0]
    yp = posp[:,1]
    zp = posp[:,2]

    xh = posh[0]
    yh = posh[1]
    zh = posh[2]

    xdel = xp - xh
    ydel = yp - yh
    zdel = zp - zh

    xp[xdel >= box/2.] = xp[xdel >= box/2.] - box
    xp[xdel < -1.*box/2.] = xp[xdel < -1. *box/2.] + box 
    yp[ydel >= box/2.] = yp[ydel >= box/2.] - box
    yp[ydel < -1.*box/2.] = yp[ydel < -1. *box/2.] + box 
    zp[zdel >= box/2.] = zp[zdel >= box/2.] - box
    zp[zdel < -1.*box/2.] = zp[zdel < -1. *box/2.] + box 
    
    posp=np.column_stack([xp,yp,zp])

    return posp

def add_ghost_particles(posc,vals,maxrad):
    #posc_ghosts = np.empty((0),float) #we aren't filling these with anything?
    #vals_ghosts = np.empty((0),float)
    posc_ghosts= posc
    vals_ghosts=vals
    print "vals before", np.shape(vals)

    x1 = -maxrad; y1 = -maxrad; z1 = -maxrad
    x2 = box + maxrad; y2 = box + maxrad; z2 = box + maxrad

    #xnew_box_size=x2-x1
    #print "xnew box size", xnew_box_size
    for i in (-1,0,1):
        for j in (-1,0,1):
            for k in (-1,0,1):
                if [i==0 and j==0 and k==0]: 
                    continue
                xp = posc[:,0] + i*box
                yp = posc[:,1] + j*box
                zp = posc[:,2] + k*box
                dm = [(xp>x1) & (xp<x2) & (yp>y1) & (yp<y2) & (zp>z1) & (zp<z2)]
                posc_new = np.column_stack([xp[dm], yp[dm],zp[dm]]); vals_new = vals[dm]
                posc_ghosts = np.concatenate((posc_ghosts,posc_new))
                vals_ghosts = np.concatenate((vals_ghosts,vals_new))
                print "i=",i,"j=",j,"k=",k

    print "vals_ghosts after ", np.shape(vals_ghosts)
    return posc_ghosts, vals_ghosts


def cull_and_center(np.ndarray posp, np.ndarray vals, np.ndarray weights, 
                    np.ndarray posh, rh,scaled_radius):

    #posp_new = periodic_bcs(posp,posh)
    #xp = posp_new[:,0]-posh[0]; yp=posp_new[:,1]-posh[1]; zp=posp_new[:,2]-posh[2]
    xp = posp[:,0]-posh[0]; yp=posp[:,1]-posh[1]; zp=posp[:,2]-posh[2]
    if (scaled_radius == True): 
        r = np.sqrt(xp**2+yp**2+zp**2)/rh
        dm = [r < search_radius]
    else:
        r = np.sqrt(xp**2+yp**2+zp**2)
        dm = [r < search_radius * rh]

    xp=xp[dm];yp=yp[dm];zp=zp[dm];vals=vals[dm];weights=weights[dm]
    posp=np.column_stack([xp,yp,zp])

    return posp,vals,weights

def precull(np.ndarray posp, np.ndarray vals, np.ndarray weights, 
            np.ndarray posh, np.ndarray rh):

    nchain = 256
    rbuff = rh.max() * search_radius

    x1 = posp[:,0].min()-1.1*rbuff; x2 = posp[:,0].max()+1.1*rbuff
    y1 = posp[:,1].min()-1.1*rbuff; y2 = posp[:,1].max()+1.1*rbuff
    z1 = posp[:,2].min()-1.1*rbuff; z2 = posp[:,2].max()+1.1*rbuff

    dx = (x2-x1) / nchain; dy = (y2-y1) / nchain; dz = (z2-z1) / nchain

    mask = np.reshape(np.zeros(nchain**3),(nchain,nchain,nchain)).astype(np.bool)
    for ih in np.arange(len(rh)):
        print ih,len(rh)
        xl = posh[ih,0] - rbuff; xh = posh[ih,0] + rbuff
        yl = posh[ih,1] - rbuff; yh = posh[ih,1] + rbuff
        zl = posh[ih,2] - rbuff; zh = posh[ih,2] + rbuff
        il = max(int((xl-x1)/dx),0); ih = min(int((xh-x1)/dx),nchain-1)
        jl = max(int((yl-y1)/dy),0); jh = min(int((yh-y1)/dy),nchain-1)
        kl = max(int((zl-z1)/dz),0); kh = min(int((zh-z1)/dz),nchain-1)
        for i in np.arange(il,ih+1):
            for j in np.arange(jl,jh+1):
                for k in np.arange(kl,kh+1):
                    mask[i,j,k] = True

    pmask = np.zeros(len(vals)).astype(np.bool)
    for ip in np.arange(len(vals)):
        if ip%10000==0: print ip,len(vals)
        x = posp[ip,0]; y = posp[ip,1]; z = posp[ip,2]
        i = int((x-x1)/dx)
        j = int((y-y1)/dy)
        k = int((z-z1)/dz)
        if mask[i,j,k]: pmask[ip] = mask[i,j,k]
                    
    posp    = posp[pmask]
    vals    = vals[pmask]
    weights = weights[pmask]

    return posp,vals,weights

def stackonhalostile(
        np.ndarray          pospi,
        np.ndarray          valsi,
        np.ndarray          poshi,
        np.ndarray            mhi,
        np.ndarray            rhi,
	np.ndarray GroupFirstSubi,
	np.ndarray           sfri,
	np.ndarray         mstari,
        it, jt, kt,ntile,volweight,mhmin, mhmax,scaled_radius):

    '''
    Parameters
	particles[nparticles][npartprops]
        halos[nphalos][nhaloprops]

    Returns
	profiles[:,nhalos]
    '''

    CHP = CompHaloProp(params.lims,params.bins)

    rpmax = rhi.max()
    rbuff=rpmax*search_radius

    xp = pospi[:,0]; yp = pospi[:,1]; zp = pospi[:,2]
    xh = poshi[:,0]; yh = poshi[:,1]; zh = poshi[:,2]

    x1=0.; x2=box; y1=0.; y2=box; z1=0.; z2=box;
    dx=(x2-x1)/ntile; dy=(y2-y1)/ntile; dz=(z2-z1)/ntile;

    x1h=it*dx; x2h=(it+1)*dx
    y1h=jt*dy; y2h=(jt+1)*dy
    z1h=kt*dz; z2h=(kt+1)*dz

    x1p=x1h-rbuff; x2p=x2h+rbuff
    y1p=y1h-rbuff; y2p=y2h+rbuff
    z1p=z1h-rbuff; z2p=z2h+rbuff

    dmp = [(xp>x1p) & (xp<x2p) & (yp>y1p) & (yp<y2p) & (zp>z1p) & (zp<z2p)]
    dmh = [(xh>x1h) & (xh<x2h) & (yh>y1h) & (yh<y2h) & (zh>z1h) & (zh<z2h) & (mhi>mhmin) & (mhi<mhmax)]

    xp=xp[dmp]; yp=yp[dmp]; zp=zp[dmp]
    xh=xh[dmh]; yh=yh[dmh]; zh=zh[dmh] 

    posp  = np.column_stack([xp,yp,zp])
    posh  = np.column_stack([xh,yh,zh])

    vals          = valsi[dmp]

    mh    	  = mhi[dmh]
    rh            = rhi[dmh] 
    GroupFirstSub = GroupFirstSubi[dmh] 
    sfr           = sfri[dmh]
    mstar         = mstari[dmh]

    pcen = np.empty((0),float)
    pval = np.empty((0),float)
    pnum = np.empty((0),float)

    nhalos=np.shape(xh)[0]
    if params.rank==0:
        print it*ntile**2+jt*ntile+kt+1,'of',ntile**3,'done, nhalos =',nhalos
    
    if nhalos == 0:
        return pcen, pval, pnum, mh, rh, nhalos, GroupFirstSub,sfr,mstar
    
    ninhalos=0
    nphalo = np.zeros(nhalos)
    
    weights = 1.0 + 0*xp

#    posp, vals, weights = precull(posp,vals,weights,posh,rh)

    for ih in np.arange(nhalos):    	
        pospc, valsc, weightsc = cull_and_center(posp,vals,weights,posh[ih],rh[ih],scaled_radius=scaled_radius)
        scale=rh[ih]
        pcenc, pvalc, pnumc = CHP.ComputeHaloProfile(pospc,valsc,weightsc,scale,volweight=volweight,scaled_radius=scaled_radius)
        pcen = np.append(pcen,pcenc)
        pval = np.append(pval,pvalc)
        pnum = np.append(pnum,pnumc)

    return pcen,pval,pnum,mh,rh,nhalos,GroupFirstSub,sfr,mstar
	
def stackonhalos(
	np.ndarray          posp,
        np.ndarray          vals,
        np.ndarray          posh,
        np.ndarray            mh,
        np.ndarray            rh,
	np.ndarray GroupFirstSub,
	np.ndarray           sfr,
	np.ndarray         mstar,
        ntile, volweight,mhmin, mhmax,scaled_radius):

    rpmax = rh.max()
    rbuff = rpmax*search_radius
    #rbuff=10000.

    posp,vals = add_ghost_particles(posp,vals,rbuff)
    
    pcen = np.empty((0),float)
    pval = np.empty((0),float)
    pnum = np.empty((0),float)
    mhpr = np.empty((0),float)
    rhpr = np.empty((0),float)
    GroupFirstSubpr=np.empty((0),float)
    sfrpr= np.empty((0),float)
    mstarpr=np.empty((0),float)
    nhalos=0
    for it in np.arange(ntile):
        for jt in np.arange(ntile):
            for kt in np.arange(ntile):

                pcenc, pvalc, pnumc, mhc, rhc, nhalosc,GroupFirstSubc,sfrc,mstarc = stackonhalostile(posp,vals,posh,mh,rh,GroupFirstSub,sfr,mstar,it,jt,kt,ntile,volweight,mhmin,mhmax,scaled_radius)   

                pcen=np.append(pcen,pcenc)
                pval=np.append(pval,pvalc)
                pnum=np.append(pnum,pnumc)
                mhpr=np.append(mhpr,  mhc)
                rhpr=np.append(rhpr,  rhc)
                GroupFirstSubpr=np.append(GroupFirstSubpr, GroupFirstSubc)
                sfrpr=np.append(sfrpr,sfrc)
                mstarpr=np.append(mstarpr,mstarc)
                nhalos += nhalosc
                
    return pcen, pval, pnum, mhpr, rhpr, nhalos, GroupFirstSubpr,sfrpr,mstarpr
		     

