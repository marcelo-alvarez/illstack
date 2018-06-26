import numpy as np

class CompHaloProp:
    def __init__(self,lims,bins,Linear=False):
        
        self.MinPos = lims[0]
        self.MaxPos = lims[1]
        self.Bins = bins

        if (Linear == True):
            self.BinSize = (self.MaxPos - self.MinPos) / self.Bins
            self.r1 = self.MinPos + self.BinSize * np.arange(self.Bins)
            self.r2 = self.MinPos + self.BinSize * (np.arange(self.Bins) + 1.0)
            self.BinCenter = self.MinPos + self.BinSize * (np.arange(self.Bins)+ 0.5)
        else:
            self.BinSize = (np.log(self.MaxPos) - np.log(self.MinPos)) / self.Bins
            self.r1 = self.MinPos * np.exp(self.BinSize * np.arange(self.Bins))
            self.r2 = self.MinPos * np.exp(self.BinSize * (np.arange(self.Bins)+ 1.0))
            self.BinCenter = self.MinPos * np.exp(self.BinSize * (np.arange(self.Bins)+ 0.5))
        
        self.radbins = np.append(self.r1,self.r2[-1])

    def ComputeHaloProfile(self,pos,quant,weight,volweight=False,stddev=False,innerbin=False):
        '''
        Returns stacked profile of a given halo
        Input: Partical position (center on Halo), Stacking quantity, Weight for average  
        Output: Bin center, Stack profile, Particle count 
        '''
        rad = np.sqrt( (pos[:,0])**2 + (pos[:,1])**2 + (pos[:,2])**2)
        Volume = 4.*np.pi/3. * (self.r2**3 - self.r1**3)

        data_qw = np.histogram(rad, bins=self.radbins, weights=quant*weight)
        data_w  = np.histogram(rad, bins=self.radbins, weights=weight)
        BinCount = np.histogram(rad, bins=self.radbins)

        if (volweight == True):
            BinValue = data_qw[0] / Volume 
        else:
            BinValue = data_qw[0] / data_w[0]

        #Add in inner bin below inner range
        if (innerbin == True):
            data_qw_inner = np.histogram(rad, bins=[0,self.radbins[0]], weights=quant*weight)
            data_w_inner  = np.histogram(rad, bins=[0,self.radbins[0]], weights=weight)
            Volume_inner = 4.*np.pi/3. * (self.r1[0]**3)

            if (volweight == True):
                BinValue[0] += data_qw_inner[0] / Volume_inner
            else:
                BinValue[0] += data_qw_inner[0] / data_w_inner[0]

        return self.BinCenter, np.nan_to_num(BinValue), BinCount[0]

    def ComputeCumulativeProfile(self,pos,quant,volweight=False,stddev=False,innerbin=True):
        '''
        Returns stacked cumulative profile of a given halo
        Input: Partical position (center on Halo), Stacking quantity
        Output: Bin center, Stack profile, Particle count
        '''
        rad = np.sqrt( (pos[:,0])**2 + (pos[:,1])**2 + (pos[:,2])**2)
        Volume = 4.*np.pi/3. * (self.r2**3 - self.r1**3)

        data_q = np.histogram(rad, bins=self.radbins, weights=quant)
        BinCount = np.histogram(rad, bins=self.radbins)

        if (volweight == True):
            BinValue = data_q[0] / Volume
        else:
            BinValue = data_q[0]

        #Add in inner bin below inner range
        if (innerbin == True):
            data_q_inner = np.histogram(rad, bins=[0,self.radbins[0]], weights=quant)
            Volume_inner = 4.*np.pi/3. * (self.r1[0]**3)
            if (volweight == True):
                BinValue[0] += data_q_inner[0] / Volume_inner
            else:
                BinValue[0] += data_q_inner[0]

        return self.BinCenter, np.nan_to_num(BinValue), BinCount[0]
