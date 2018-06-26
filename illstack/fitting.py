import numpy as np
import emcee
import time

class fitting_mcmc():
    def __init__(self,nwalkers=100,chains=5000,burnin=50,verbose=False):
        self.nwalkers = nwalkers
        self.chains = chains
        self.burnin = burnin

    def run_mcmc(self,Ndim,P0,lnprob,args,verbose=False):
        '''
        the arg format follows the function lnprob
        e.g. lnprob(theat,a,b,c) then arg = (a,b,c) 
        '''
        pos = [P0 + P0*1e-1*np.random.randn(Ndim) for i in range(self.nwalkers)]

        start = time.clock()

        sampler = emcee.EnsembleSampler(nwalkers,Ndim,lnprob, args = args)
        sampler.run_mcmc(pos,self.chains)
        elapsed1 = (time.clock() - start)
        if (verbose == True):
            print elapsed1

        samples = sampler.chain[:,self.burnin:,:].reshape((-1,Ndim))    
        ans_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                           zip(*np.percentile(samples, [16, 50, 84],axis=0)))

        return ans_mcmc,samples
