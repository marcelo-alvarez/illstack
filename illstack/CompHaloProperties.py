import numpy as np

class CompHaloProp:
    def __init__(self,lims,bins,Linear=False):
        
        self.MinPos = lims[0]
        self.MaxPos = lims[1]
        self.Bins = bins

        if (Linear == True):
            print ("Linear binning")
            self.BinSize = (self.MaxPos - self.MinPos) / self.Bins
            self.r1 = self.MinPos + self.BinSize * np.arange(self.Bins)
            self.r2 = self.MinPos + self.BinSize * (np.arange(self.Bins) + 1.0)
            self.BinCenter = self.MinPos + self.BinSize * (np.arange(self.Bins)+ 0.5)
        else:
            print ("Log binning")
            self.BinSize = (np.log(self.MaxPos) - np.log(self.MinPos)) / self.Bins
            self.r1 = self.MinPos * np.exp(self.BinSize * np.arange(self.Bins))
            self.r2 = self.MinPos * np.exp(self.BinSize * (np.arange(self.Bins)+ 1.0))
            self.BinCenter = self.MinPos * np.exp(self.BinSize * (np.arange(self.Bins)+ 0.5))
        
        self.radbins = np.append(self.r1,self.r2[-1])

    def ComputeHaloProfile(self,pos,quant,weight,volweight=False,StdDev=False):

        rad = np.sqrt( (pos[0,:])**2 + (pos[1,:])**2 + (pos[2,:])**2)
        Volume = 4.*np.pi/3. * (self.r2**3 - self.r1**3)

        data_qw = np.histogram(rad, bins=self.radbins, weights=quant*weight)
        data_w  = np.histogram(rad, bins=self.radbins, weights=weight)
        data_count = np.histogram(rad, bins=self.radbins)

        if (volweight == True):
            BinValue = data_qw[0] / Volume
        else:
            BinValue = data_qw[0] / data_w[0]

        return self.BinCenter, BinValue, BinCount
